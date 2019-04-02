import argparse
from basics import *
from gnomad_hail import *
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gnomAD_liftover")
logger.setLevel(logging.INFO)


def checkpoint(ht: hl.Table, dtype: str) -> hl.Table:
    """
    Creates a "checkpoint" by writing Table out and reading back in
    :param Table ht: Table to be evaluated and saved
    :param str dtype: Data type (exomes or genomes)
    :return: Table after writing out and reading back in
    :rtype: Table
    """

    if dtype == 'exome':
        ht.write('gs://seqr-datasets/methods_dev/test_data/kc/gnomAD_exome_sites.ht', overwrite=True)
        ht = hl.read_table('gs://seqr-datasets/methods_dev/test_data/kc/gnomAD_exome_sites.ht')

    else:
        ht.write('gs://seqr-datasets/methods_dev/test_data/kc/gnomAD_genome_sites.ht', overwrite=True)
        ht = hl.read_table('gs://seqr-datasets/methods_dev/test_data/kc/gnomAD_genome_sites.ht')

    return ht


def lift_ht(ht: hl.Table, dtype: str, rg: hl.genetics.ReferenceGenome) -> hl.Table:
    """
    Lifts gnomAD release table from b37 to b38
    :param Table ht: gnomAD release Table
    :param str dtype: Data type (exome or genome)
    :param ReferenceGenome rg: Reference genome, b38
    :return: gnomAD release Table lifted to b38
    :rtype: Table
    """

    # get reference genomes and load liftover chain file
    #rg37 = hl.get_reference('GRCh37')
    #rg38 = rg38 = hl.get_reference('GRCh38')
    #rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

    # annotate release file with liftover coordinates
    ht = ht.annotate(new_locus = hl.liftover(ht.locus, rg, include_strand=True),
                     old_locus=ht.locus
                    )
    ht = ht.key_by(locus = ht.new_locus.result)
    ht = ht.annotate(no_target = hl.is_defined(ht.locus))
    ht = checkpoint(ht, dtype)
    return ht


def reverse_complement(allele: str):
    """
    Returns the reverse complement of an allele
    :param str allele: Allele to be flipped
    :return: Reverse complement of input allele
    :rtype: str
    """
 
    return (hl.switch(allele)
            .when('A', 'T')
            .when('T', 'A')
            .when('G', 'C')
            .when('C', 'G')
            .default(allele))


def snp_mismatch(ht: hl.Table, dtype: str, rg: hl.genetics.ReferenceGenome) -> hl.Table:
    """
    Annotates mismatches between reference allele and allele in reference fasta
    :param Table ht: Table of SNPs to be annotated
    :param str dtype: Data type (exome or genome)
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
                                (reverse_complement(ht.alleles[0]) != hl.get_sequence(ht.locus.contig, ht.locus.position, reference_genome=rg)),
                                (ht.alleles[0] != hl.get_sequence(ht.locus.contig, ht.locus.position, reference_genome=rg))
                                )
                )

    #ht = checkpoint(ht, dtype)
    return ht
 
def check_mismatch(ht: hl.Table) -> dict:
    """
    Checks for mismatches between reference allele and allele in reference fasta
    :param Table ht: Table to be checked
    :return: Dictionary containing mismatch type ('total' or 'minus' as key) and umber of mismatches as value
    :rtype: dict
    """

    mismatch = {}
    ht = ht.filter(ht.reference_mismatch)
    mismatch['total'] = ht.count()
    mismatch['minus'] =  ht.filter(ht.new_locus.is_negative_strand).count()
    return mismatch


def main(args):

    hl.init()
    
    # get reference genomes and load liftover chain file
    rg37 = hl.get_reference('GRCh37')
    rg38 = rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38)

    # add fasta sequence to rg38
    rg38.add_sequence(
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz', 
            'gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai'
            )
    
    logger.info('Reading in release ht')
    if args.type == 'exome':
        dtype = 'exome'
        ht = hl.read_table(public_exomes_ht_path(split=True, version=CURRENT_RELEASE))
    else:
        dtype = 'genome'
        ht = hl.read_table(public_genomes_ht_path(split=True, version=CURRENT_RELEASE))

    logger.info('Working on gnomAD {} release ht'.format(dtype))
    #logger.info('Filtering to chr21 for testing')
    #ht = ht.filter(ht.locus.contig == '21')

    logger.info('Lifting ht to b38')
    ht = lift_ht(ht, dtype, rg38)

    logger.info('Checking SNPs for reference mismatches')
    ht = snp_mismatch(ht, dtype, rg38)
    mismatch = check_mismatch(ht)
    logger.info('{} reference mismatches in SNPs'.format(mismatch['total']))
    logger.info('{} mismatches on minus strand'.format(mismatch['minus']))


if __name__ == '__main__':

    # Create argument parser
    parser = argparse.ArgumentParser(description='This script lifts gnomad ht to b38')
    parser.add_argument('-t', '--type', help='Data type', choices=['exome', 'genome'], required=True)
    #parser.add_argument('-o', '--overwrite', help='Overwrite pre-existing data', action='store_true', default=True)
    args = parser.parse_args()
    main(args)
