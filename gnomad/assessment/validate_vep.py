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


def get_variant_count_in_protein_coding_gene(
    release: str("v4.0"), dataset, vep_version: str("105"), test=False
):
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

    startTime = datetime.now()

    # import interval information of protein-coding genes
    if test == True:
        i = hl.import_table(
            "gs://gnomad-qin/test_3pc_genes.tsv", delimiter="\t", impute=True
        )
    else:
        i = hl.import_table(
            f"gs://gnomad-qin/ensembl_v{vep_version}_pc_genes.tsv",
            delimiter="\t",
            impute=True,
        )

    # define the path of the VEP-annotated HT
    if release == "v4.0" and dataset == "exomes":
        t = hl.read_table(
            "gs://gnomad/v4.0/annotations/exomes/gnomad.exomes.v4.0.vep.ht"
        )
    elif release == "v4.0" and dataset == "genomes":
        t = hl.read_table(
            "gs://gnomad/v4.0/annotations/genomes/gnomad.genomes.v4.0.vep.ht"
        )
    elif release == "v3.1" and dataset == "genomes":
        t = hl.read_table(
            "gs://gnomad/release/3.1.4/ht/genomes/gnomad.genomes.v3.1.4.sites.ht"
        )
    else:
        logger.warning(
            "currently only support v4.0 exomes, v4.0 genomes and v3.1 genomes"
        )

    # convert the interval information to locus_interval in Hail format
    i = i.key_by(
        interval=hl.locus_interval(
            hl.literal("chr") + hl.str(i.chr), i.start, i.end, reference_genome="GRCh38"
        )
    )

    i = i.checkpoint(
        f"gs://gnomad-tmp-4day/interval_info_{release}_{dataset}_{vep_version}.ht",
        overwrite=True,
    )

    # join the interval information with the VEP-annotated HT by locus
    i = i.annotate(interval_copy=i.interval)
    t = t.annotate(the_interval_or_NA=i[t.locus].interval_copy)  # join by locus
    t = t.filter(hl.is_defined(t.the_interval_or_NA))
    t = t.annotate(the_interval=t.the_interval_or_NA)

    # count the number of total variants and "protein-coding" variants in each interval
    re = t.group_by(t.the_interval).aggregate(
        all_variants=hl.agg.count(),
        variants_in_pcg=hl.agg.count_where(
            t.vep.transcript_consequences.biotype.contains("protein_coding")
        ),
    )

    re = re.checkpoint(
        f"gs://gnomad-tmp-4day/interval_count_{release}_{dataset}_{vep_version}.ht",
        overwrite=True,
    )

    i = i.annotate(
        all_variants=re[i.interval].all_variants,
        variants_in_pcg=re[i.interval].variants_in_pcg,
    )

    i = i.checkpoint(
        f"gs://gnomad-tmp-4day/interval_info_with_count_{release}_{dataset}_{vep_version}.ht",
        overwrite=True,
    )

    # output the number of genes and the name of genes
    na_genes = i.filter(hl.is_missing(i.all_variants)).gene_stable_ID.collect()
    logger.info(
        f"{len(na_genes)} gene(s) have no overlapped intervals with this gnomAD release"
    )

    variant0_genes = i.filter(i.all_variants == 0).gene_stable_ID.collect()
    logger.info(
        f"{len(variant0_genes)} gene(s) have their intervals overlapped but have no"
        " variants found in their defined intervals"
    )

    non_pcg_genes = i.filter(
        (i.all_variants != 0) & (i.variants_in_pcg == 0)
    ).gene_stable_ID.collect()
    logger.info(
        f"{len(non_pcg_genes)} gene(s) have variants found in overlapped interval"
        " but no variants annotated as protein-coding"
        f" biotype in their intervals, including {non_pcg_genes}"
    )

    partial_pcg_genes = i.filter(
        (i.all_variants != 0)
        & (i.variants_in_pcg != 0)
        & (i.all_variants != i.variants_in_pcg)
    ).gene_stable_ID.collect()
    logger.info(
        f"{len(partial_pcg_genes)} gene(s) have a subset of variants annotated as"
        " protein-coding biotype in their defined intervals"
    )

    logger.info(f"Time elapsed: {datetime.now() - startTime}")


def main(args):
    """Validate number of variants by VEP in protein-coding genes."""
    get_variant_count_in_protein_coding_gene(
        args.release, args.dataset, args.vep_version, args.test
    )


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
