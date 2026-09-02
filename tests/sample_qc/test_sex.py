"""Tests for the sex ploidy adjustment functions in gnomad.sample_qc.sex."""

import hail as hl
import pytest

from gnomad.sample_qc.sex import adjust_sex_ploidy, adjusted_sex_ploidy_expr


class TestAdjustSexPloidy:
    """Test that `adjust_sex_ploidy` matches `adjusted_sex_ploidy_expr`."""

    @pytest.fixture
    def mt(self):
        """
        Build an MT covering every branch of the sex ploidy adjustment.

        Loci: an autosome, chrX PAR1 / non-PAR (x2) / PAR2, chrY PAR1 / non-PAR.
        Samples: XX, XY, lowercase "xy" (matched case-insensitively), an
        unrecognized karyotype and a missing one. Genotypes cycle through
        hom-ref, het, hom-var and missing.
        """
        loci = [
            hl.locus("chr1", 1000, reference_genome="GRCh38"),
            hl.locus("chrX", 20000, reference_genome="GRCh38"),
            hl.locus("chrX", 5000000, reference_genome="GRCh38"),
            hl.locus("chrX", 100000000, reference_genome="GRCh38"),
            hl.locus("chrY", 20000, reference_genome="GRCh38"),
            hl.locus("chrY", 5000000, reference_genome="GRCh38"),
            hl.locus("chrX", 155800000, reference_genome="GRCh38"),
        ]
        samples = [
            ("s1", "XX"),
            ("s2", "XY"),
            ("s3", "xy"),
            ("s4", "XXY"),
            ("s5", None),
        ]
        gts = [(0, 0), (0, 1), (1, 1), None]
        entries = []
        for l_idx, locus in enumerate(loci):
            for s_idx, (s, _) in enumerate(samples):
                gt = gts[(l_idx + s_idx) % len(gts)]
                entries.append(
                    {
                        "locus": locus,
                        "s": s,
                        "GT": hl.call(*gt) if gt is not None else hl.missing(hl.tcall),
                    }
                )
        mt = hl.Table.parallelize(
            entries, hl.tstruct(locus=hl.tlocus("GRCh38"), s=hl.tstr, GT=hl.tcall)
        ).to_matrix_table(row_key=["locus"], col_key=["s"])
        karyotype = hl.literal(dict(samples), hl.tdict(hl.tstr, hl.tstr))
        return mt.annotate_cols(sex_karyotype=karyotype.get(mt.s))

    @staticmethod
    def _entries(mt, gt_field):
        return sorted(
            (e.locus.contig, e.locus.position, e.s, e[gt_field])
            for e in mt.entries().collect()
        )

    def test_matches_expr_version(self, mt):
        """The in-place adjustment matches `adjusted_sex_ploidy_expr` entry for entry."""
        expected = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype)
        )
        actual = adjust_sex_ploidy(mt, mt.sex_karyotype)
        assert self._entries(actual, "GT") == self._entries(expected, "GT")
        # The temporary flag fields must not leak into the output.
        assert set(actual.row) == set(mt.row)
        assert set(actual.col) == set(mt.col)

    def test_adjusts_expected_entries(self, mt):
        """Spot-check every branch of the adjustment on known entries."""
        actual = self._entries(adjust_sex_ploidy(mt, mt.sex_karyotype), "GT")
        by_key = {(c, p, s): gt for c, p, s, gt in actual}
        # Genotype of sample s_idx at locus l_idx is gts[(l_idx + s_idx) % 4].
        # Autosome and PAR genotypes are untouched (XY, PAR1 and PAR2).
        assert by_key[("chr1", 1000, "s2")] == hl.Call([0, 1])
        assert by_key[("chrX", 20000, "s2")] == hl.Call([1, 1])
        assert by_key[("chrX", 155800000, "s3")] == hl.Call([0, 0])
        # XY on non-PAR X/Y: het -> missing, hom -> haploid, missing stays missing.
        assert by_key[("chrX", 100000000, "s3")] is None  # het, lowercase "xy"
        assert by_key[("chrX", 100000000, "s2")] == hl.Call([0])
        assert by_key[("chrX", 5000000, "s2")] is None
        assert by_key[("chrY", 5000000, "s2")] == hl.Call([1])
        # XX on Y (PAR or not) -> missing regardless of the genotype.
        assert by_key[("chrY", 20000, "s1")] is None
        assert by_key[("chrY", 5000000, "s1")] is None
        # Unrecognized or missing karyotypes are left alone.
        assert by_key[("chrX", 5000000, "s4")] == hl.Call([0, 1])
        assert by_key[("chrY", 5000000, "s5")] == hl.Call([0, 1])

    def test_gt_field(self, mt):
        """`gt_field` selects a non-default genotype entry field."""
        mt = mt.rename({"GT": "LGT"})
        expected = mt.annotate_entries(
            LGT=adjusted_sex_ploidy_expr(mt.locus, mt.LGT, mt.sex_karyotype)
        )
        actual = adjust_sex_ploidy(mt, mt.sex_karyotype, gt_field="LGT")
        assert self._entries(actual, "LGT") == self._entries(expected, "LGT")
