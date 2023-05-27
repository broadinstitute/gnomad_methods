"""script to validate VEP annotation on protein-coding genes."""
import argparse
import logging
from datetime import datetime

import hail as hl

from gnomad.resources.grch38.reference_data import ensembl_interval

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def count_variant_per_interval(
    ht, release: str("v4.0"), dataset: str("exomes"), vep_version: str("105")
) -> hl.Table:
    """
    Count how many variants annotated by VEP in each protein-coding gene according to its interval information.

    :param ht: hl.Table(), ensembl interval HT
    :param release: str, "v4.0" or "v3.1"
    :param dataset: str, "exomes" or "genomes"
    :param vep_version: str, "101" or "105"
    :return: vep_ht: hl.Table(), VEP-annotated HT, only selected the vep.transcript_consequences field
    """
    logger.info(
        "importing vep_HT on dataset: gnomAD_%s_%s_vep%s", release, dataset, vep_version
    )
    # check the input parameters
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
    # TODO: changed this part and relevant function to import the vep HT
    # resources, instead of hard-coding the path
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

    count_ht = count_ht.checkpoint("gs://gnomad-tmp-4day/count_tmp.ht", overwrite=True)

    ht = ht.annotate(**count_ht[ht.gene_stable_ID])

    logger.info("checkpointing the count by interval HT: ")
    ht.checkpoint(
        f"gs://gnomad-tmp-4day/validate_vep/count_interval_{release}_{dataset}_{vep_version}.ht",
        overwrite=True,
    )
    # TODO: change the path after PR review

    logger.info("Reporting genes without variants annotated: ")
    na_genes = ht.filter(
        hl.is_missing(ht.variants_in_pcg) | (ht.variants_in_pcg == 0)
    ).gene_stable_ID.collect()

    logger.info(
        "%s gene(s) have no variants annotated as protein-coding in Biotype. It is"
        " likely these genes are not covered by this gnomAD release. These genes"
        " are: %s",
        len(na_genes),
        na_genes,
    )

    partial_pcg_genes = ht.filter(
        (ht.all_variants != 0)
        & (ht.variants_in_pcg != 0)
        & (ht.all_variants != ht.variants_in_pcg)
    ).gene_stable_ID.collect()
    logger.info(
        "%s gene(s) have a subset of variants annotated as protein-coding biotype in"
        " their defined intervals",
        {len(partial_pcg_genes)},
    )

    return ht


def main(args):
    """Validate number of variants by VEP in protein-coding genes."""
    hl.init(log="/tmp/hail.log", tmp_dir="gs://gnomad-tmp-4day")
    startTime = datetime.now()

    logger.info("importing and parsing interval file: ")
    ht = ensembl_interval.ht()
    # TODO: ask Julia to copy HT and raw *.tsv to resources location

    logger.info("counting variants in the interval of each protein-coding gene: ")

    count_variant_per_interval(ht, args.release, args.dataset, args.vep_version)

    logger.info("Time elapsed: %s", datetime.now() - startTime)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "--release",
        default="v4.0",
        help="gnomAD release version",
        required=True,
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
