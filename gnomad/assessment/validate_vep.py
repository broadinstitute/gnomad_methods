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


def count_variant_per_interval(ht, vep_ht) -> hl.Table:
    """
    Count total number of variants and variants annotated to be protein-coding in each interval.

    :param ht: hl.Table(), containing interval information of protein-coding genes
    :param vep_ht: hl.Table(), VEP-annotated HT, only selected the vep.transcript_consequences field
    :return: hl.Table(), interval information with variant count
    """
    # join two tables by interval and locus as index key
    vep_ht = vep_ht.annotate(
        interval_annotations=ht.index(vep_ht.locus, all_matches=True)
    )

    vep_ht = vep_ht.filter(hl.is_defined(vep_ht.interval_annotations))

    vep_ht = vep_ht.annotate(
        gene_stable_ID=vep_ht.interval_annotations.gene_stable_ID,
        biotype=vep_ht.transcript_consequences.biotype,
    )

    # select only the gene_stable_ID and biotype to save space
    vep_ht = vep_ht.select("gene_stable_ID", "biotype")

    # explode the vep_ht by gene_stable_ID
    vep_ht = vep_ht.explode(vep_ht.gene_stable_ID)

    # count the number of total variants and "protein-coding" variants in each interval
    count_ht = vep_ht.group_by(vep_ht.gene_stable_ID).aggregate(
        all_variants=hl.agg.count(),
        variants_in_pcg=hl.agg.count_where(vep_ht.biotype.contains("protein_coding")),
    )

    count_ht.checkpoint("gs://gnomad-tmp-4day/count_tmp.ht", overwrite=True)

    ht = ht.annotate(**count_ht[ht.gene_stable_ID])

    na_genes = ht.filter(hl.is_missing(ht.variants_in_pcg) | (ht.variants_in_pcg == 0))
    logger.info(
        f"{len(na_genes.gene_stable_ID.collect())} gene(s) have no variants annotated"
        " as protein-coding in Biotype. It is likely these genes are not covered by"
        f" gnomAD data. These genes are: {na_genes.gene_stable_ID.collect()}"
    )

    partial_pcg_genes = ht.filter(
        (ht.all_variants != 0)
        & (ht.variants_in_pcg != 0)
        & (ht.all_variants != ht.variants_in_pcg)
    ).gene_stable_ID.collect()
    logger.info(
        f"{len(partial_pcg_genes)} gene(s) have a subset of variants annotated as"
        " protein-coding biotype in their defined intervals"
    )

    return ht


def import_filter_vep_data(
    release: str("v4.0"), dataset: str("exomes"), vep_version: str("105")
) -> hl.Table:
    """
    Count how many variants annotated by VEP in each protein-coding gene according to its interval information.

    :param release: str, "v4.0" or "v3.1"
    :param dataset: str, "exomes" or "genomes"
    :param vep_version: str, "101" or "105"
    :return: vep_ht: hl.Table(), VEP-annotated HT, only selected the vep.transcript_consequences field
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

    return vep_ht


def import_parse_interval(interval_file) -> hl.Table:
    """
    Import and parse interval of protein-coding genes to HT.

    Downloaded from Ensembl Archive for 101 & 105.
    :param interval_file: str, the name of the interval file
    """
    ht = hl.import_table(
        f"gs://gnomad-qin/{interval_file}.tsv",
        delimiter="\t",
        min_partitions=100,
        impute=True,
    )

    ht = ht.key_by(
        interval=hl.locus_interval(
            hl.literal("chr") + hl.str(ht.chr),
            ht.start,
            ht.end,
            reference_genome="GRCh38",
        )
    )
    ht = ht.checkpoint(
        f"gs://gnomad-tmp-4day/qin/{interval_file}.ht", _read_if_exists=True
    )
    return ht


def main(args):
    """Validate number of variants by VEP in protein-coding genes."""
    hl.init(log="/tmp/hail.log", tmp_dir="gs://gnomad-tmp-4day")
    startTime = datetime.now()

    logger.info("importing and parsing interval file: {}".format(args.interval_file))
    ht = import_parse_interval(args.interval_file)

    logger.info(
        "importing vep_HT on dataset:"
        f" gnomAD_{args.release}_{args.dataset}_vep{args.vep_version}"
    )
    vep_ht = import_filter_vep_data(args.release, args.dataset, args.vep_version)

    logger.info(f"counting variants in the interval of each protein-coding gene: ")
    ht = count_variant_per_interval(ht, vep_ht)

    logger.info("writing out the count by interval HT: ")
    ht.checkpoint(
        f"gs://gnomad-tmp-4day/count_interval_{args.release}_{args.dataset}_{args.vep_version}.ht",
        overwrite=True,
    )

    logger.info(f"Time elapsed: {datetime.now() - startTime}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--release", default="v4.0", help="gnomAD release version", required=True
    )
    parser.add_argument(
        "--interval-file", help="prefix of interval file", required=True
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

    main(parser.parse_args())
