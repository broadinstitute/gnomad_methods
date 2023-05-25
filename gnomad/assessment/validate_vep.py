"""script to validate VEP annotation on protein-coding genes."""
import argparse
import logging
from datetime import datetime

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def count_variant_per_interval(i_ht: hl.Table(), vep_ht: hl.Table()) -> hl.Table():
    """
    Count how many variants in each interval.

    :param i_ht: hl.Table(), interval information of protein-coding genes, downloaded from Ensembl Archive for 101 & 105
    :param vep_ht: hl.Table(), VEP-annotated HT, only selected the vep.transcript_consequences field
    :return: hl.Table(), interval information with variant count
    """
    # convert the interval information to locus_interval in Hail format
    i_ht = i_ht.key_by(
        interval=hl.locus_interval(
            hl.literal("chr") + hl.str(i_ht.chr),
            i_ht.start,
            i_ht.end,
            reference_genome="GRCh38",
        )
    )

    # join the interval information with the VEP-annotated HT by locus
    i_ht = i_ht.annotate(interval_copy=i_ht.interval)
    vep_ht = vep_ht.annotate(
        the_interval_or_NA=i_ht[vep_ht.locus].interval_copy
    )  # join by locus
    vep_ht = vep_ht.filter(hl.is_defined(vep_ht.the_interval_or_NA))
    vep_ht = vep_ht.annotate(the_interval=vep_ht.the_interval_or_NA)

    # count the number of total variants and "protein-coding" variants in each interval
    re = vep_ht.group_by(vep_ht.the_interval).aggregate(
        all_variants=hl.agg.count(),
        variants_in_pcg=hl.agg.count_where(
            vep_ht.transcript_consequences.biotype.contains("protein_coding")
        ),
    )
    i_ht = i_ht.annotate(
        all_variants=re[i_ht.interval].all_variants,
        variants_in_pcg=re[i_ht.interval].variants_in_pcg,
    )
    na_genes = i_ht.filter(hl.is_missing(i_ht.all_variants))
    logger.info(
        f"{len(na_genes).gene_stable_ID.collect()} gene(s) have no variants annotated"
        " as protein-coding in biotype."
    )
    return i_ht


def import_filter_vep_data(
    release: str("v4.0"), dataset, vep_version: str("105"), test=False
) -> hl.Table():
    """
    Count how many variants annotated by VEP in each protein-coding gene according to its interval information.

    :param release: str, "v4.0" or "v3.1"
    :param dataset: str, "exomes" or "genomes"
    :param vep_version: str, "101" or "105"
    :param test: bool, if True, only run on 3 protein-coding genes
    :return:
    """
    if release == "v4.0" and vep_version == "101" and dataset == "exomes":
        logger.warning(
            "gnomAD v4.0 EXOMES does not have v101 VEP annotation, please use 105"
            " instead"
        )

    if release == "v4.0" and vep_version == "101" and dataset == "genomes":
        logger.warning(
            "gnomAD v4.0 GENOMES does not have v101 VEP annotation, please use 105"
            " instead"
        )

    if release == "v3.1" and vep_version == "105" and dataset == "genomes":
        logger.warning(
            "gnomAD v3.1 GENOMES does not have v105 VEP annotation, please use 101"
            " instead"
        )

    if release not in ["v4.0", "v3.1"]:
        logger.warning("release should be v4.0 or v3.1 only")

    if vep_version not in ["101", "105"]:
        logger.warning("vep_version should be 101 or 105 only")

    # import interval information of protein-coding genes
    if test == True:
        i_ht = hl.import_table(
            "gs://gnomad-qin/test_3pc_genes.tsv", delimiter="\t", impute=True
        )
    else:
        i_ht = hl.import_table(
            f"gs://gnomad-qin/ensembl_v{vep_version}_pc_genes.tsv",
            delimiter="\t",
            impute=True,
        )

    # define the path of the VEP-annotated HT
    if release == "v4.0" and dataset == "exomes":
        vep_ht = hl.read_table(
            "gs://gnomad/v4.0/annotations/exomes/gnomad.exomes.v4.0.vep.ht"
        )
    elif release == "v4.0" and dataset == "genomes":
        vep_ht = hl.read_table(
            "gs://gnomad/v4.0/annotations/genomes/gnomad.genomes.v4.0.vep.ht"
        )
    elif release == "v3.1" and dataset == "genomes":
        vep_ht = hl.read_table(
            "gs://gnomad/release/3.1.4/ht/genomes/gnomad.genomes.v3.1.4.sites.ht"
        )
    else:
        logger.warning(
            "currently only support v4.0 exomes, v4.0 genomes and v3.1 genomes"
        )

    vep_ht = vep_ht.select(vep_ht.vep.transcript_consequences)

    return i_ht, vep_ht


def filter_na_genes(i_ht: hl.Table()) -> hl.Table():
    """Filter genes with no variants annotated as protein-coding in biotype.

    :param i_ht: hl.Table, interval HT with count of variants per gene.
    """
    na_genes = i_ht.filter(hl.is_missing(i_ht.all_variants))
    logger.info(
        "Filtering genes with no variants annotated as protein-coding in biotype."
    )
    return na_genes


def main(args):
    """Validate number of variants by VEP in protein-coding genes."""
    hl.init(log="/tmp/hail.log", tmp_dir="gs://gnomad-tmp-4day")
    startTime = datetime.now()
    logger.info(
        "importing vep_HT on dataset:"
        f" gnomAD_{args.release}_{args.dataset}_vep{args.vep_version}"
    )
    i_ht, vep_ht = import_filter_vep_data(
        args.release, args.dataset, args.vep_version, args.test
    )
    logger.info(f"counting variants in each protein-coding gene: 1st round")
    i_ht = count_variant_per_interval(i_ht, vep_ht)
    logger.info(
        f"filtering genes with no variants annotated as protein-coding in biotype"
    )
    na_genes = filter_na_genes(i_ht)
    logger.info(f"counting variants in each protein-coding gene: 2nd round")
    na_genes_counted = count_variant_per_interval(na_genes)
    logger.info(f"combining 2 rounds of counting")
    i_ht = i_ht.annotate(
        all_variants=na_genes_counted[i_ht.interval].all_variants,
        variants_in_pcg=na_genes_counted[i_ht.interval].variants_in_pcg,
    )
    i_ht.checkpoint(
        f"gs://gnomad-qin/validate_vep_results/interval_info_with_count_{args.release}_{args.dataset}_vep{args.vep_version}.ht",
        overwrite=True,
    )
    logger.info(f"Time elapsed: {datetime.now() - startTime}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--release", default="v4.0", help="gnomAD release version", required=True
    )
    parser.add_argument(
        "--dataset",
        default="exomes",
        help="gnomAD dataset, either exomes or genomes",
        required=True,
    )
    parser.add_argument(
        "--vep-version",
        default="105",
        help="VEP version, either 101 or 105",
        required=True,
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="whether to run on test interval file, default is false",
    )
    main(parser.parse_args())
