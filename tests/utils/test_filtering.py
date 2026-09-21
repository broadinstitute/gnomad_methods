"""Tests for the filtering utility module."""

import hail as hl
import pytest

from gnomad.utils.filtering import filter_to_autosomes


class TestFilterToAutosomes:
    """Test the filter_to_autosomes function."""

    @pytest.mark.parametrize(
        "build, loci",
        [
            ("GRCh38", ["chr1:100", "chr22:100", "chrX:100", "chrY:100", "chrM:100"]),
            ("GRCh37", ["1:100", "22:100", "X:100", "Y:100", "MT:100"]),
        ],
    )
    def test_filter_to_autosomes(self, build: str, loci: list) -> None:
        """Test that sex and MT contigs are removed and autosomes retained."""
        ht = hl.Table.parallelize(
            [{"locus": hl.parse_locus(l, reference_genome=build)} for l in loci],
            hl.tstruct(locus=hl.tlocus(build)),
            key="locus",
        )
        ht = filter_to_autosomes(ht)
        assert [str(r.locus) for r in ht.collect()] == loci[:2]

    def test_matrixtable(self) -> None:
        """Test that a MatrixTable input is filtered on rows."""
        ht = hl.Table.parallelize(
            [
                {"locus": hl.parse_locus(l, reference_genome="GRCh38")}
                for l in ["chr2:5", "chrX:5"]
            ],
            hl.tstruct(locus=hl.tlocus("GRCh38")),
            key="locus",
        )
        mt = hl.MatrixTable.from_rows_table(ht)
        mt = filter_to_autosomes(mt)
        assert mt.count_rows() == 1
