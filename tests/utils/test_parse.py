"""Tests for the parse utility module."""

import pytest
import hail as hl
from gnomad.utils.parse import parse_variant

class TestParseVariant:
    """Test the parse_variant function."""

    @pytest.mark.parametrize(
        "variant_str, contig, position, ref, alt, build, expected",
        [
            ("1-1000-A-T", None, None, None, None, None, hl.Struct(locus=hl.Locus("1", 1000, "GRCh37"), alleles=["A", "T"])),
            ("chr1-1000-A-T", None, None, None, None, None, hl.Struct(locus=hl.Locus("chr1", 1000, "GRCh38"), alleles=["A", "T"])),
            (None, "1", 1000, "A", "T", None, hl.Struct(locus=hl.Locus("1", 1000, "GRCh37"), alleles=["A", "T"])),
            (None, "chr1", 1000, "A", "T", None, hl.Struct(locus=hl.Locus("chr1", 1000, "GRCh38"), alleles=["A", "T"])),
            (None, "1", 1000, "A", "T", "GRCh37", hl.Struct(locus=hl.Locus("1", 1000, "GRCh37"), alleles=["A", "T"])),
            ("1:1000:A:T", None, None, None, None, None, hl.Struct(locus=hl.Locus("1", 1000, "GRCh37"), alleles=["A", "T"])),
        ],
    )
    def test_parse_variant(
        self,
        variant_str: str,
        contig: str,
        position: int,
        ref: str,
        alt: str,
        build: str,
        expected: hl.expr.StructExpression,
    ) -> None:
        """
        Test valid parameters for the `parse_variant` function.

        :param variant_str: Variant string.
        :param contig: Chromosome of the variant.
        :param position: Variant position.
        :param ref: Reference allele.
        :param alt: Alternate allele.
        :param build: Reference genome build.
        :param expected: Expected result.
        :return: None.
        """
        result = hl.eval(parse_variant(variant_str, contig, position, ref, alt, build))
        assert result == expected

    @pytest.mark.parametrize(
        "variant_str, contig, position, ref, alt, build",
        [
            (None, None, None, None, None, None),
            ("invalid_variant", None, None, None, None, None),
            (None, "1", None, "A", "T", None),
        ],
    )
    def test_parse_variant_invalid(
        self,
        variant_str: str,
        contig: str,
        position: int,
        ref: str,
        alt: str,
        build: str,
    ) -> None:
        """
        Test invalid parameters for the `parse_variant` function.

        :param variant_str: Variant string.
        :param contig: Chromosome of the variant.
        :param position: Variant position.
        :param ref: Reference allele.
        :param alt: Alternate allele.
        :param build: Reference genome build.
        :return: None.
        """
        with pytest.raises(ValueError):
            parse_variant(variant_str, contig, position, ref, alt, build)
