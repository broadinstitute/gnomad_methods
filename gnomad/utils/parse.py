"""This module contains utility functions for general parsing."""

import logging
from typing import List, Optional, Union

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("parse_utils")
logger.setLevel(logging.INFO)


def parse_variant(
    variant_str: Optional[str] = None,
    contig: Optional[str] = None,
    position: Optional[int] = None,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
    build: Optional[str] = None,
) -> hl.expr.StructExpression:
    """
    Create a Struct with the locus and alleles from a variant string or contig, position, ref, and alt.

    :param variant_str: Variant string in the format contig-position-ref-alt or
        contig:position:ref:alt.
    :param contig: Chromosome of the variant.
    :param position: Variant position.
    :param ref: Reference allele.
    :param alt: Alternate allele.
    :param build: Reference genome build. If not provided, will infer from the variant
        string or contig. If 'chr' is present in the contig, will assume GRCh38,
        otherwise GRCh37.
    :return: Struct with the locus and alleles.
    """
    if not variant_str and not all([contig, position, ref, alt]):
        raise ValueError(
            "Either `variant_str` must be provided or all of `contig`, `position`, "
            "`ref`, and `alt`."
        )

    if not build:
        build = "GRCh37"
        if (variant_str and variant_str.startswith("chr")) or (
            contig and contig.startswith("chr")
        ):
            build = "GRCh38"

        logger.info("No build provided. Assuming build: %s", build)

    try:
        if variant_str and ":" not in variant_str:
            contig, position, ref, alt = variant_str.split("-")
        if all([contig, position, ref, alt]):
            variant_str = f"{contig}:{position}:{ref}:{alt}"

        return hl.parse_variant(variant_str, reference_genome=build)

    except BaseException:
        raise ValueError(
            f"Invalid variant format: {variant_str}. Valid formats: \n"
            f"  contig-position-ref-alt \n"
            f"  contig:position:ref:alt"
        )


def parse_locus_intervals(
    intervals: Union[str, List[str]],
    reference_genome: Optional[str] = None,
) -> Union[hl.expr.IntervalExpression, List[hl.expr.IntervalExpression]]:
    """
    Parse interval(s) string into Hail locus intervals.

    .. note::

        If no reference genome is provided, the function will infer the reference genome
        build from the interval string. If the interval string starts with 'chr', it will
        assume GRCh38, otherwise GRCh37.
    :param intervals: Interval string or list of interval strings. The interval string
        format has to be "contig:start-end", e.g.,"1:1000-2000" (GRCh37) or
        "chr1:1000-2000" (GRCh38).
    :return: Hail locus interval(s).
    """
    is_str = isinstance(intervals, str)
    intervals = [intervals] if is_str else intervals

    if not reference_genome:
        reference_genome = "GRCh38" if intervals[0].startswith("chr") else "GRCh37"
        logger.info(
            "No reference genome provided. Assuming reference genome: %s",
            reference_genome,
        )

    if reference_genome == "GRCh38" and any(
        [not i.startswith("chr") for i in intervals]
    ):
        raise ValueError("Interval must start with 'chr' for GRCh38 reference genome.")

    intervals = [
        hl.parse_locus_interval(i, reference_genome=reference_genome) for i in intervals
    ]

    return intervals[0] if is_str else intervals
