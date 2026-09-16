"""Tests for the constraint utility module."""

from typing import List

import hail as hl
import pytest

from gnomad.utils.constraint import assemble_constraint_context_ht


def _context_ht(build: str, loci: List[str]) -> hl.Table:
    """Build a minimal split context Table with heptamer context and VEP annotation."""
    rows = [
        {
            "locus": hl.parse_locus(l, reference_genome=build),
            "alleles": ["C", "T"],
            "context": "AAACGAA",
            "was_split": False,
            "vep": hl.struct(
                most_severe_consequence="missense_variant",
                transcript_consequences=[
                    hl.struct(
                        transcript_id="ENST1",
                        gene_id="ENSG1",
                        consequence_terms=["missense_variant"],
                    )
                ],
            ),
        }
        for l in loci
    ]
    return hl.Table.parallelize(
        rows,
        hl.tstruct(
            locus=hl.tlocus(build),
            alleles=hl.tarray(hl.tstr),
            context=hl.tstr,
            was_split=hl.tbool,
            vep=hl.tstruct(
                most_severe_consequence=hl.tstr,
                transcript_consequences=hl.tarray(
                    hl.tstruct(
                        transcript_id=hl.tstr,
                        gene_id=hl.tstr,
                        consequence_terms=hl.tarray(hl.tstr),
                    )
                ),
            ),
        ),
        key=["locus", "alleles"],
    )


class TestAssembleConstraintContextHt:
    """Test the assemble_constraint_context_ht function."""

    @pytest.mark.parametrize(
        "build, keep, drop",
        [
            (
                "GRCh38",
                ["chr1:100", "chr22:100", "chrX:100", "chrY:100"],
                ["chrM:100", "chr1_KI270706v1_random:100"],
            ),
            ("GRCh37", ["1:100", "22:100", "X:100", "Y:100"], ["MT:100"]),
        ],
    )
    def test_filters_to_primary_contigs(
        self, build: str, keep: List[str], drop: List[str]
    ) -> None:
        """Test that only autosomes and sex contigs are retained."""
        ht = assemble_constraint_context_ht(_context_ht(build, keep + drop))
        assert sorted(str(r.locus) for r in ht.collect()) == sorted(keep)

    def test_annotations(self) -> None:
        """Test context trimming, strand collapse, and mutation type annotations."""
        ht = assemble_constraint_context_ht(_context_ht("GRCh38", ["chr1:100"]))
        row = ht.collect()[0]
        assert row.context == "ACG"
        assert row.ref == "C" and row.alt == "T"
        assert not row.was_flipped
        assert row.cpg and row.mutation_type == "CpG"
        assert row.vep.transcript_consequences[0].most_severe_consequence == (
            "missense_variant"
        )
        assert "consequence_terms" not in row.vep.transcript_consequences[0]
