import argparse
from gnomad_hail.resources import *
from gnomad_hail.utils.generic import flip_base 
import logging
from os.path import dirname, basename
import sys


logging.basicConfig(format="%(asctime)s (%(name)s %(lineno)s): %(message)s", datefmt='%m/%d/%Y %I:%M:%S %p')
logger = logging.getLogger("liftover")
logger.setLevel(logging.INFO)


def get_checkpoint_path(gnomad: bool, data_type: str, path: str) -> str:
    """
    Creates a checkpoint path for Table
    :param bool gnomad: Whether data is gnomAD data
    :param str data_type: Data type (exomes or genomes for gnomAD; not used otherwise)
    :param str path: Path to input Table/MatrixTable (if data is not gnomAD data)
    :return: Output checkpoint path
    :rtype: str
    """
    
    if gnomad:
        return f'gs://gnomad/liftover/release/2.1.1/ht/{data_type}/gnomad.{data_type}.r2.1.1.liftover.ht'
    else:
        out_name = basename(path).split('.')[0]
        return f'{dirname(path)}/{out_name}_lift.ht'


def lift_data(t: Union[hl.MatrixTable, hl.Table], gnomad: bool, data_type: str, path: str, rg: hl.genetics.ReferenceGenome) -> Union[hl.MatrixTable, hl.Table]:
    """
    Lifts input Table or MatrixTable from one reference build to another

    :param Table/MatrixTable t: Table or MatrixTable
    :param bool gnomad: Whether data is gnomAD data
    :param str data_type: Data type (exomes or genomes for gnomAD; not used otherwise)
    :param str path: Path to input Table/MatrixTable (if data is not gnomAD data)
    :param ReferenceGenome rg: Reference genome
    :return: Table or MatrixTablewith liftover annotations
    :rtype: Table or MatrixTable
    """

    logger.info('Annotate input with liftover coordinates')
    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(new_locus=hl.liftover(t.locus, rg, include_strand=True),
                            old_locus=ht.locus)
        num_no_target = t.aggregate_rows(hl.agg.count_where(hl.is_missing(t.new_locus)))
        logger.info(f'Filtering out {num_no_target} sites that failed to liftover')
        t = t.filter_rows(hl.is_defined(t.new_locus)) 
        t = t.key_rows_by(locus=t.new_locus.result, alleles=t.alleles)

    else:
        t = t.annotate(new_locus=hl.liftover(t.locus, rg, include_strand=True),
                     old_locus=ht.locus
                    )
        num_no_target = t.aggregate(hl.agg.count_where(hl.is_missing(t.new_locus)))
        logger.info(f'Filtering out {num_no_target} sites that failed to liftover')
        t = t.filter(hl.is_defined(t.new_locus)) 
        t = t.key_by(locus=t.new_locus.result, alleles=t.alleles)

    logger.info('Writing out lifted over data')
    t = t.checkpoint(get_checkpoint_path(gnomad, data_type, path), overwrite=True)
    return ht


def annotate_snp_mismatch(t: Union[hl.MatrixTable, hl.Table], data_type: str, rg: hl.genetics.ReferenceGenome) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotates mismatches between reference allele and allele in reference fasta

    Assumes input Table/MatrixTable has t.new_locus annotation
    :param Table/MatrixTable t: Table/MatrixTable of SNPs to be annotated
    :param str data_type: Data type (exomes or genomes for gnomAD; not used otherwise)
    :param ReferenceGenome rg: Reference genome with fasta sequence loaded
    :return: Table annotated with mismatches between reference allele and allele in fasta
    :rtype: Table
    """
 
    logger.info('Filtering to SNPs')
    if isinstance(t, hl.MatrixTable):
        t = t.filter_rows(hl.is_snp(t.alleles[0], t.alleles[1]))

        # check if reference allele matches what is in reference fasta
        # for snps on negative strand, make sure reverse complement of ref allele matches what is in ref fasta
        t = t.annotate_rows(
            reference_mismatch=hl.cond(
                                t.new_locus.is_negative_strand,
                                (flip_base(ht.alleles[0]) != hl.get_sequence(t.locus.contig, t.locus.position, reference_genome=rg)),
                                (t.alleles[0] != hl.get_sequence(t.locus.contig, t.locus.position, reference_genome=rg))
                                ))
    else:
        t = t.filter(hl.is_snp(t.alleles[0], t.alleles[1]))
        t = t.annotate(
            reference_mismatch=hl.cond(
                                t.new_locus.is_negative_strand,
                                (flip_base(ht.alleles[0]) != hl.get_sequence(t.locus.contig, t.locus.position, reference_genome=rg)),
                                (t.alleles[0] != hl.get_sequence(t.locus.contig, t.locus.position, reference_genome=rg))
                                ))
    return t

 
def check_mismatch(t: Union[hl.MatrixTable, hl.Table]) -> hl.expr.expressions.StructExpression:
    """
    Checks for mismatches between reference allele and allele in reference fasta
    :param Table ht: Table to be checked
    :return: StructExpression containing counts for mismatches and count for all variants on negative strand
    :rtype: StructExpression
    """

    mismatch = ht.aggregate(hl.struct(total_variants=hl.agg.count(),
                                total_mismatch=hl.agg.count_where(ht.reference_mismatch),
                                negative_strand=hl.agg.count_where(ht.new_locus.is_negative_strand),
                                negative_strand_mismatch=hl.agg.count_where(ht.new_locus.is_negative_strand & ht.reference_mismatch)
                                )
                            )
    return mismatch


def main(args):

    hl.init('/liftover.log')
    
    if args.gnomad:
        gnomad = True
        path = None

        if args.exomes:
            data_type = 'exomes'
        if args.genomes:
            data_type = 'genomes'

        logger.info('Working on gnomAD {} release ht'.format(data_type))
        logger.info('Reading in release ht')
        t = get_gnomad_public_data(data_type, split=True, version=CURRENT_RELEASE)
        logger.info('Variants in release ht: {}'.format(ht.count()))

    else:
        data_type = None
        gnomad = False
   
        if args.ht:
            path = args.ht
            t = hl.read_table(args.ht)
        if args.mt:
            path = args.mt
            t = hl.read_table(args.mt)
   
    logger.info('Checking if input data has been split') 
    rows = t.row
    if 'was_split' not in rows:
        if isinstance(t, hl.MatrixTable):
            t = split_multi_hts(t)
        else:
            t = split_multi(t)

    logger.info('Inferring build of input')
    build = hl.get_reference_genome(t.locus).name
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

    if args.test: 
        logger.info('Filtering to chr21 for testing')
        if build == 38:
            contig = 'chr21'
        else:
            contig = '21'
        if isinstance(t, hl.MatrixTable):
            t = t.filter_rows(t.locus.contig == contig)
        else:
            t = t.filter(t.locus.contig == contig)

    logger.info(f'Lifting data to {target.name}')
    ht = lift_ht(ht, gnomad, data_type, path, target)
        
    logger.info('Checking SNPs for reference mismatches')
    ht = annotate_snp_mismatch(ht, data_type, target)
    
    if isinstance(t, hl.MatrixTable):
        mismatch = check_mismatch(t.rows())
    else:
        mismatch = check_mismatch(t)
    logger.info('{} total SNPs'.format(mismatch['total_variants']))
    logger.info('{} SNPs on minus strand'.format(mismatch['negative_strand']))
    logger.info('{} reference mismatches in SNPs'.format(mismatch['total_mismatch']))
    logger.info('{} mismatches on minus strand'.format(mismatch['negative_strand_mismatch']))


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='This script lifts a ht from one build to another')
    parser.add_argument('--mt', help='Full path to MatrixTable to liftover. Specify only if not using --gnomad flag')
    parser.add_argument('--ht', help='Full path to Table to liftover. Specify only if not using --gnomad flag.')
    parser.add_argument('-g', '--gnomad', help='Liftover table is one of the gnomAD releases', action='store_true')
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
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if not args.gnomad and (int(args.mt) + int(args.ht) != 1):
        sys.exit('Error: One and only one of --mt or --ht must be specified')

    main(args)
