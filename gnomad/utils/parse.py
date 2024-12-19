"""This module contains utility functions for general parsing."""
import logging
from typing import Optional

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
        if (variant_str and variant_str.startswith("chr")) or (contig and contig.startswith("chr")):
            build = "GRCh38"

        logger.info("No build provided. Assuming build: %s", build)

    try:
        if variant_str and ":" not in variant_str:
            contig, position, ref, alt = variant_str.split("-")
        if all([contig, position, ref, alt]):
            variant_str = f"{contig}:{position}:{ref}:{alt}"

        return hl.parse_variant(variant_str, reference_genome=build)

    except:
        raise ValueError(
            f"Invalid variant format: {variant_str}. Valid formats: \n"
            f"  contig-position-ref-alt \n"
            f"  contig:position:ref:alt"
        )
