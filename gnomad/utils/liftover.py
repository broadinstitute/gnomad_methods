import argparse
from typing import Union
import hail as hl
from gnomad.resources.grch37.gnomad import public_release
from gnomad.utils.generic import flip_base, get_reference_genome
import logging
from os.path import dirname, basename
import sys


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("liftover")
logger.setLevel(logging.INFO)


def get_checkpoint_path(gnomad: bool, data_type: str, path: str, is_table: bool) -> str:
    """
    Creates a checkpoint path for Table

    :param gnomad: Whether data is gnomAD data
    :param data_type: Data type (exomes or genomes for gnomAD; not used otherwise)
    :param path: Path to input Table/MatrixTable (if data is not gnomAD data)
    :param is_table: Whether data is a Table
    :return: Output checkpoint path
    """
    
    if gnomad:
        return f'gs://gnomad/liftover/release/2.1.1/ht/{data_type}/gnomad.{data_type}.r2.1.1.liftover.ht'
    else:
        out_name = basename(path).split('.')[0]
        out = f'{dirname(path)}/{out_name}_lift'
        return f'{out}.ht'if is_table else f'{out}.mt'


def get_liftover_genome(t: Union[hl.MatrixTable, hl.Table]) -> list:
    """
    Infers genome build of input data and assumes destination build. Prepares to liftover to destination genome build

    :param t: Input Table or MatrixTable
    :return: List of source build (with liftover chain added) and destination build (with sequence loaded)
    """

    logger.info('Inferring build of input')
    build = get_reference_genome(t.locus).name

    logger.info('Loading reference genomes, adding chain file, and loading fasta sequence for destination build')
    if build == 'GRCh38':
        source = hl.get_reference('GRCh38')
        target = hl.get_reference('GRCh37')
        chain = 'gs://hail-common/references/grch38_to_grch37.over.chain.gz'
        target.add_sequence(
            'gs://hail-common/references/human_g1k_v37.fasta.gz',
            'gs://hail-common/references/human_g1k_v37.fasta.fai'
            ) 
    else:
        source = hl.get_reference('GRCh37')
        target = hl.get_reference('GRCh38')
        chain = 'gs://hail-common/references/grch37_to_grch38.over.chain.gz'
        target.add_sequence(
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz', 
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai'
            )

    source.add_liftover(chain, target)
    return [source, target]


def lift_data(t: Union[hl.MatrixTable, hl.Table], gnomad: bool, data_type: str, path: str, rg: hl.genetics.ReferenceGenome,
                overwrite: bool) -> Union[hl.MatrixTable, hl.Table]:
    """
    Lifts input Table or MatrixTable from one reference build to another

    :param t: Table or MatrixTable
    :param gnomad: Whether data is gnomAD data
    :param data_type: Data type (exomes or genomes for gnomAD; not used otherwise)
    :param path: Path to input Table/MatrixTable (if data is not gnomAD data)
    :param rg: Reference genome
    :param overwrite: Whether to overwrite data
    :return: Table or MatrixTablewith liftover annotations
    """

    logger.info('Annotating input with liftover coordinates')
    liftover_expr = {
        'new_locus': hl.liftover(t.locus, rg, include_strand=True),
        'old_locus': t.locus
    }
    t = t.annotate(**liftover_expr) if isinstance(t, hl.Table) else t.annotate_rows(**liftover_expr)

    no_target_expr = hl.agg.count_where(hl.is_missing(t.new_locus))
    num_no_target = t.aggregate(no_target_expr) if isinstance(t, hl.Table) else t.aggregate_rows(no_target_expr)
    
    logger.info(f'Filtering out {num_no_target} sites that failed to liftover')
    keep_expr = hl.is_defined(t.new_locus)
    t = t.filter(keep_expr) if isinstance(t, hl.Table) else t.filter_rows(keep_expr)
    
    row_key_expr = {
        'locus': t.new_locus.result,
        'alleles': t.alleles
    }
    t = t.key_by(**row_key_expr) if isinstance(t, hl.Table) else t.key_rows_by(**row_key_expr)

    logger.info('Writing out lifted over data')
    t = t.checkpoint(get_checkpoint_path(gnomad, data_type, path, isinstance(t, hl.Table)), overwrite=overwrite)
    return t


def annotate_snp_mismatch(t: Union[hl.MatrixTable, hl.Table], data_type: str, rg: hl.genetics.ReferenceGenome) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotates mismatches between reference allele and allele in reference fasta

    Assumes input Table/MatrixTable has t.new_locus annotation

    :param t: Table/MatrixTable of SNPs to be annotated
    :param data_type: Data type (exomes or genomes for gnomAD; not used otherwise)
    :param rg: Reference genome with fasta sequence loaded
    :return: Table annotated with mismatches between reference allele and allele in fasta
    """
 
    logger.info('Filtering to SNPs')
    snp_expr = hl.is_snp(t.alleles[0], t.alleles[1])
    t = t.filter(snp_expr) if isinstance(t, hl.Table) else t.filter_rows(snp_expr)

    mismatch_expr = {
        'reference_mismatch': hl.cond(t.new_locus.is_negative_strand, 
                                    (flip_base(t.alleles[0]) != hl.get_sequence(t.locus.contig, t.locus.position, reference_genome=rg)),
                                    (t.alleles[0] != hl.get_sequence(t.locus.contig, t.locus.position, reference_genome=rg)))
    }
    logger.info('Checking if reference allele matches what is in reference fasta')
    logger.info('For SNPs on the negative strand, make sure the reverse complement of the ref alleles matches what is in the ref fasta')
    return t.annotate(**mismatch_expr) if isinstance(t, hl.Table) else t.annotate_rows(**mismatch_expr)
    
 
def check_mismatch(ht: hl.Table) -> hl.expr.expressions.StructExpression:
    """
    Checks for mismatches between reference allele and allele in reference fasta

    :param ht: Table to be checked
    :return: StructExpression containing counts for mismatches and count for all variants on negative strand
    """

    mismatch = ht.aggregate(hl.struct(total_variants=hl.agg.count(),
                                total_mismatch=hl.agg.count_where(ht.reference_mismatch),
                                negative_strand=hl.agg.count_where(ht.new_locus.is_negative_strand),
                                negative_strand_mismatch=hl.agg.count_where(ht.new_locus.is_negative_strand & ht.reference_mismatch)
                                )
                            )
    return mismatch


def main(args):

    hl.init(log='/liftover.log')
    
    if args.gnomad:
        gnomad = True
        path = None
        
        if args.exomes:
            data_type = 'exomes'
        if args.genomes:
            data_type = 'genomes'

        logger.info('Working on gnomAD {} release ht'.format(data_type))
        logger.info('Reading in release ht')
        t = public_release(data_type).ht()
        logger.info('Variants in release ht: {}'.format(t.count()))

    else:
        data_type = None
        gnomad = False
   
        if args.ht:
            path = args.ht
            t = hl.read_table(args.ht)
        if args.mt:
            path = args.mt
            t = hl.read_matrix_table(args.mt)
   
    logger.info('Checking if input data has been split') 
    if 'was_split' not in t.row:
        t = hl.split_multi(t) if isinstance(t, hl.Table) else hl.split_multi_hts(t) 

    logger.info('Preparing reference genomes for liftover')
    source, target = get_liftover_genome(t)
    
    if args.test: 
        logger.info('Filtering to chr21 for testing')
        if source.name == 'GRCh38':
            contig = 'chr21'
        else:
            contig = '21'
        t = hl.filter_intervals(t, [hl.parse_locus_interval(contig, reference_genome=source.name)])

    logger.info(f'Lifting data to {target.name}')
    t = lift_data(t, gnomad, data_type, path, target, args.overwrite)
        
    logger.info('Checking SNPs for reference mismatches')
    t = annotate_snp_mismatch(t, data_type, target)
   
    mismatch = check_mismatch(t) if isinstance(t, hl.Table) else check_mismatch(t.rows()) 
    logger.info('{} total SNPs'.format(mismatch['total_variants']))
    logger.info('{} SNPs on minus strand'.format(mismatch['negative_strand']))
    logger.info('{} reference mismatches in SNPs'.format(mismatch['total_mismatch']))
    logger.info('{} mismatches on minus strand'.format(mismatch['negative_strand_mismatch']))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='This script lifts a ht from one build to another')
    parser.add_argument('--mt', help='Full path to MatrixTable to liftover. Specify only if not using --gnomad flag')
    parser.add_argument('--ht', help='Full path to Table to liftover. Specify only if not using --gnomad flag.')
    parser.add_argument('-g', '--gnomad', help='Liftover table is one of the gnomAD releases', action='store_true')
    parser.add_argument('-o', '--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument(
            '--exomes', 
            help='Data type is exomes. One of --exomes or --genomes is required if --gnomad is specified.', 
            action='store_true'
            )
    parser.add_argument(
            '--genomes', 
            help='Data type is genomes. One of --exomes or --genomes is required if --gnomad is specified.', 
            action='store_true'
            )
    parser.add_argument('-t', '--test', help='Filter to chr21 (for code testing purposes)', action='store_true')
    args = parser.parse_args()

    if args.gnomad and (int(args.exomes) + int(args.genomes) != 1):
        sys.exit('Error: One and only one of --exomes or --genomes must be specified with --gnomad flag')

    if args.mt and args.ht:
        sys.exit('Error: One and only one of --mt or --ht must be specified')

    main(args)
