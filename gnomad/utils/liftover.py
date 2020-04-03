import logging
from os.path import basename, dirname
from typing import Union

import hail as hl
from gnomad.utils.reference_genome import get_reference_genome

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
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
        return f"gs://gnomad/liftover/release/2.1.1/ht/{data_type}/gnomad.{data_type}.r2.1.1.liftover.ht"
    else:
        out_name = basename(path).split(".")[0]
        out = f"{dirname(path)}/{out_name}_lift"
        return f"{out}.ht" if is_table else f"{out}.mt"


def get_liftover_genome(t: Union[hl.MatrixTable, hl.Table]) -> list:
    """
    Infers genome build of input data and assumes destination build. Prepares to liftover to destination genome build

    :param t: Input Table or MatrixTable
    :return: List of source build (with liftover chain added) and destination build (with sequence loaded)
    """

    logger.info("Inferring build of input")
    build = get_reference_genome(t.locus).name

    logger.info(
        "Loading reference genomes, adding chain file, and loading fasta sequence for destination build"
    )
    if build == "GRCh38":
        source = hl.get_reference("GRCh38")
        target = hl.get_reference("GRCh37")
        chain = "gs://hail-common/references/grch38_to_grch37.over.chain.gz"
        target.add_sequence(
            "gs://hail-common/references/human_g1k_v37.fasta.gz",
            "gs://hail-common/references/human_g1k_v37.fasta.fai",
        )
    else:
        source = hl.get_reference("GRCh37")
        target = hl.get_reference("GRCh38")
        chain = "gs://hail-common/references/grch37_to_grch38.over.chain.gz"
        target.add_sequence(
            "gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz",
            "gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai",
        )

    source.add_liftover(chain, target)
    return [source, target]


def lift_data(
    t: Union[hl.MatrixTable, hl.Table],
    gnomad: bool,
    data_type: str,
    path: str,
    rg: hl.genetics.ReferenceGenome,
    overwrite: bool,
) -> Union[hl.MatrixTable, hl.Table]:
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

    logger.info("Annotating input with liftover coordinates")
    liftover_expr = {
        "new_locus": hl.liftover(t.locus, rg, include_strand=True),
        "old_locus": t.locus,
    }
    t = (
        t.annotate(**liftover_expr)
        if isinstance(t, hl.Table)
        else t.annotate_rows(**liftover_expr)
    )

    no_target_expr = hl.agg.count_where(hl.is_missing(t.new_locus))
    num_no_target = (
        t.aggregate(no_target_expr)
        if isinstance(t, hl.Table)
        else t.aggregate_rows(no_target_expr)
    )

    logger.info(f"Filtering out {num_no_target} sites that failed to liftover")
    keep_expr = hl.is_defined(t.new_locus)
    t = t.filter(keep_expr) if isinstance(t, hl.Table) else t.filter_rows(keep_expr)

    row_key_expr = {"locus": t.new_locus.result, "alleles": t.alleles}
    t = (
        t.key_by(**row_key_expr)
        if isinstance(t, hl.Table)
        else t.key_rows_by(**row_key_expr)
    )

    logger.info("Writing out lifted over data")
    t = t.checkpoint(
        get_checkpoint_path(gnomad, data_type, path, isinstance(t, hl.Table)),
        overwrite=overwrite,
    )
    return t


def annotate_snp_mismatch(
    t: Union[hl.MatrixTable, hl.Table], rg: hl.genetics.ReferenceGenome
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Annotates mismatches between reference allele and allele in reference fasta

    Assumes input Table/MatrixTable has t.new_locus annotation

    :param t: Table/MatrixTable of SNPs to be annotated
    :param rg: Reference genome with fasta sequence loaded
    :return: Table annotated with mismatches between reference allele and allele in fasta
    """

    logger.info("Filtering to SNPs")
    snp_expr = hl.is_snp(t.alleles[0], t.alleles[1])
    t = t.filter(snp_expr) if isinstance(t, hl.Table) else t.filter_rows(snp_expr)

    mismatch_expr = {
        "reference_mismatch": hl.cond(
            t.new_locus.is_negative_strand,
            (
                hl.reverse_complement(t.alleles[0])
                != hl.get_sequence(
                    t.locus.contig, t.locus.position, reference_genome=rg
                )
            ),
            (
                t.alleles[0]
                != hl.get_sequence(
                    t.locus.contig, t.locus.position, reference_genome=rg
                )
            ),
        )
    }
    logger.info("Checking if reference allele matches what is in reference fasta")
    logger.info(
        "For SNPs on the negative strand, make sure the reverse complement of the ref alleles matches what is in the ref fasta"
    )
    return (
        t.annotate(**mismatch_expr)
        if isinstance(t, hl.Table)
        else t.annotate_rows(**mismatch_expr)
    )


def check_mismatch(ht: hl.Table) -> hl.expr.expressions.StructExpression:
    """
    Checks for mismatches between reference allele and allele in reference fasta

    :param ht: Table to be checked
    :return: StructExpression containing counts for mismatches and count for all variants on negative strand
    """

    mismatch = ht.aggregate(
        hl.struct(
            total_variants=hl.agg.count(),
            total_mismatch=hl.agg.count_where(ht.reference_mismatch),
            negative_strand=hl.agg.count_where(ht.new_locus.is_negative_strand),
            negative_strand_mismatch=hl.agg.count_where(
                ht.new_locus.is_negative_strand & ht.reference_mismatch
            ),
        )
    )
    return mismatch


def liftover_using_gnomad_map(ht, data_type):
    """
    Liftover a gnomAD table using already-established liftover file. Warning: shuffles!

    :param ht: Input Hail table
    :param data_type: one of "exomes" or "genomes" which to map across
    :return: Lifted over table
    """
    from gnomad.resources.grch37.gnomad import liftover

    lift_ht = liftover(data_type).ht()
    ht = ht.key_by(original_locus=ht.locus, original_alleles=ht.alleles).drop(
        "locus", "alleles"
    )
    return lift_ht.annotate(
        **ht[(lift_ht.original_locus, lift_ht.original_alleles)]
    ).key_by("locus", "alleles")
