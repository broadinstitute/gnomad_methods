# noqa: D100

import logging
from typing import Tuple, Union

import hail as hl

from gnomad.utils.reference_genome import add_reference_sequence, get_reference_genome

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


GRCH37_to_GRCH38_CHAIN = "gs://hail-common/references/grch37_to_grch38.over.chain.gz"
"""
Path to chain file required to lift data from GRCh37 to GRCh38.
"""

GRCH38_TO_GRCH37_CHAIN = "gs://hail-common/references/grch38_to_grch37.over.chain.gz"
"""
Path to chain file required to lift data from GRCh38 to GRCh37.
"""


def get_liftover_genome(
    t: Union[hl.MatrixTable, hl.Table]
) -> Tuple[hl.genetics.ReferenceGenome, hl.genetics.ReferenceGenome]:
    """
    Infer reference genome build of input data and assume destination reference genome build.

    Adds liftover chain to source reference genome and sequence to destination reference genome.
    Returns tuple containing both reference genomes in preparation for liftover.

    :param t: Input Table or MatrixTable.
    :return: Tuple of source reference genome (with liftover chain added)
        and destination reference genome (with sequence loaded)
    """
    logger.info("Inferring reference genome of input...")
    input_build = get_reference_genome(t.locus).name
    source = hl.get_reference(input_build)

    logger.info("Loading fasta sequence for destination build...")
    if input_build == "GRCh38":
        target = hl.get_reference("GRCh37")
        chain = GRCH38_TO_GRCH37_CHAIN

    else:
        target = hl.get_reference("GRCh38")
        chain = GRCH37_to_GRCH38_CHAIN

    logger.info("Adding liftover chain to input build...")
    if source.has_liftover(target):
        logger.warning(
            "Source reference build %s already has a chain file: %s! Using whichever"
            " chain has already been added.",
            source.name,
            source._liftovers,
        )
    else:
        source.add_liftover(chain, target)

    return (source, add_reference_sequence(target))


def liftover_expr(
    locus: hl.expr.LocusExpression,
    alleles: hl.expr.ArrayExpression,
    destination_reference: hl.ReferenceGenome,
) -> hl.expr.StructExpression:
    """
    Generate struct liftover expression.

    Struct contains:
        - locus: Liftover coordinates
        - alleles: Liftover alleles
        - original_locus: Locus prior to liftover
        - original_alleles: Alleles prior to liftover
        - locus_fail_liftover: Whether the locus failed liftover
        - ref_allele_mismatch: Whether the allele at index 0 of alleles (lifted over reference allele)
            doesn't match the allele at that position in the destination reference

    :param locus: Input locus.
    :param alleles: Input alleles.
    :param destination_reference: Desired reference genome build for liftover.
    :return: Struct containing expressions for lifted over locus/alleles as well as original locus/alleles.
    """
    lifted_over_locus = hl.liftover(locus, destination_reference, include_strand=True)
    lifted_over_alleles = alleles.map(
        lambda a: hl.if_else(
            lifted_over_locus.is_negative_strand, hl.reverse_complement(a), a
        )
    )

    return hl.struct(
        new_locus=lifted_over_locus.result,
        new_alleles=lifted_over_alleles,
        original_locus=locus,
        original_alleles=alleles,
        locus_fail_liftover=hl.is_missing(lifted_over_locus),
        ref_allele_mismatch=lifted_over_locus.result.sequence_context()
        != lifted_over_alleles[0],
    )


def default_lift_data(
    t: Union[hl.MatrixTable, hl.Table],
    remove_failed_sites: bool = True,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Lift input Table or MatrixTable from one reference build to another.

    :param t: Table or MatrixTable.
    :return: Table or MatrixTable with liftover annotations.
    """
    logger.info("Inferring input reference and destination reference...")
    _, target_genome = get_liftover_genome(t)

    logger.info("Annotating input data with liftover coordinates...")
    t = (
        t.annotate(**liftover_expr(t.locus, t.alleles, target_genome))
        if isinstance(t, hl.Table)
        else t.annotate_rows(**liftover_expr(t.locus, t.alleles, target_genome))
    )

    no_target_expr = hl.agg.count_where(t.locus_fail_liftover)
    num_no_target = (
        t.aggregate(no_target_expr)
        if isinstance(t, hl.Table)
        else t.aggregate_rows(no_target_expr)
    )

    if remove_failed_sites:
        logger.info("Filtering out %d sites that failed to liftover...", num_no_target)
        keep_expr = ~t.locus_fail_liftover
        t = t.filter(keep_expr) if isinstance(t, hl.Table) else t.filter_rows(keep_expr)

    row_key_expr = {"locus": t.new_locus, "alleles": t.new_alleles}
    return (
        t.key_by(**row_key_expr)
        if isinstance(t, hl.Table)
        else t.key_rows_by(**row_key_expr)
    )


def liftover_using_gnomad_map(ht: hl.Table, data_type: str):
    """
    Liftover a gnomAD v2 table using already-established liftover file.

    .. note::
        This function shuffles!

    :param ht: Input Hail Table.
    :param data_type: Which gnomAD data type to map across. One of "exomes" or "genomes".
    :return: Lifted over Table
    """
    from gnomad.resources.grch37.gnomad import liftover

    logger.warning("This function will trigger a shuffle! Pre-emptibles may not work.")
    lift_ht = liftover(data_type).ht()
    ht = ht.key_by(original_locus=ht.locus, original_alleles=ht.alleles).drop(
        "locus", "alleles"
    )
    return lift_ht.annotate(
        **ht[(lift_ht.original_locus, lift_ht.original_alleles)]
    ).key_by("locus", "alleles")
