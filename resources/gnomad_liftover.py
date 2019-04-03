import argparse
from gnomad_hail.resources import *
from gnomad_hail.utils import *
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gnomAD_liftover")
logger.setLevel(logging.INFO)


def get_checkpoint_path(data_type: str) -> str:
    """
    Creates a checkpoint path for Table
    :param str data_type: Data type (exomes or genomes)
    :return: Output checkpoint path
    :rtype: str
    """

    return f'gs://seqr-datasets/methods_dev/test_data/kc/gnomAD_{data_type}_sites.ht'


def lift_ht(ht: hl.Table, data_type: str, rg: hl.genetics.ReferenceGenome) -> hl.Table:
    """
    Lifts gnomAD release table from b37 to b38
    :param Table ht: gnomAD release Table
    :param str data_type: Data type (exome or genome)
    :param ReferenceGenome rg: Reference genome, b38
    :return: gnomAD release Table lifted to b38
    :rtype: Table
    """

    # annotate release file with liftover coordinates
    ht = ht.annotate(new_locus = hl.liftover(ht.locus, rg, include_strand=True),
                     old_locus=ht.locus
                    )
    ht = ht.key_by(locus = ht.new_locus.result, alleles = ht.alleles)
    ht = ht.checkpoint(get_checkpoint_path(data_type), overwrite=True)
    return ht


def snp_mismatch(ht: hl.Table, data_type: str, rg: hl.genetics.ReferenceGenome) -> hl.Table:
    """
    Annotates mismatches between reference allele and allele in reference fasta
    :param Table ht: Table of SNPs to be annotated
    :param str data_type: Data type (exome or genome)
    :param ReferenceGenome rg: GRCh38 reference genome with fasta sequence loaded
    :return: Table annotated with mismatches between reference allele and allele in fasta
    :rtype: Table
    """
 
    # filter to snps
    ht = ht.filter(hl.is_snp(ht.alleles[0], ht.alleles[1]))
    snp_count = ht.count()
    logger.info('{} SNPs in table'.format(snp_count))

    # check if reference allele matches what is in reference fasta
    # for snps on negative strand, make sure reverse complement of ref allele matches what is in reference fasta
    ht = ht.annotate(
        reference_mismatch = hl.cond(
                                ht.new_locus.is_negative_strand,
                                (generic.flip_base(ht.alleles[0]) != hl.get_sequence(ht.locus.contig, ht.locus.position, reference_genome=rg)),
                                (ht.alleles[0] != hl.get_sequence(ht.locus.contig, ht.locus.position, reference_genome=rg))
                                )
                )

    return ht

 
def check_mismatch(ht: hl.Table) -> hl.expr.expressions.StructExpression:
    """
    Checks for mismatches between reference allele and allele in reference fasta
    :param Table ht: Table to be checked
    :return: StructExpression containing counts for mismatches and count for all variants on negative strand
    :rtype: StructExpression
    """

    mismatch = ht.aggregate(hl.struct(total_mismatch=hl.agg.count_where(ht.reference_mismatch),
                                negative_strand=hl.agg.count_where(ht.new_locus.is_negative_strand),
                                negative_strand_mismatch=hl.agg.count_where(ht.new_locus.is_negative_strand & ht.reference_mismatch)
                                )
                            )
    return mismatch


def main(args):

    hl.init()
    
    # get reference genomes and load liftover chain file
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

    # add fasta sequence to rg38
    rg38.add_sequence(
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz', 
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai'
            )
    
    if args.exomes:
        data_type = 'exomes'
    if args.genomes:
        data_type = 'genomes'

    logger.info('Working on gnomAD {} release ht'.format(data_type))
    logger.info('Reading in release ht')
    ht = get_gnomad_public_data(data_type, split=True, version=CURRENT_RELEASE)
    
    #logger.info('Filtering to chr21 for testing')
    #ht = ht.filter(ht.locus.contig == '21')

    logger.info('Lifting ht to b38')
    ht = lift_ht(ht, data_type, rg38)

    ht = hl.read_table(f'gs://seqr-datasets/methods_dev/test_data/kc/gnomAD_{data_type}_sites.ht')
    
    logger.info('Checking SNPs for reference mismatches')
    ht = snp_mismatch(ht, data_type, rg38)
    mismatch = check_mismatch(ht)
    logger.info('{} reference mismatches in SNPs'.format(mismatch['total_mismatch']))
    logger.info('{} mismatches on minus strand'.format(mismatch['negative_strand_mismatch']))


if __name__ == '__main__':

    # Create argument parser
    parser = argparse.ArgumentParser(description='This script lifts gnomAD release ht to b38')
    parser.add_argument('--exomes', help='Data type is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Data type is genomes. One of --exomes or --genomes is required.', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    main(args)
