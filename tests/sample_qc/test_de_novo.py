"""Test suite for de novo mutation functions."""

import hail as hl
import pytest

from gnomad.sample_qc.relatedness import (
    calculate_de_novo_post_prob,
    default_get_de_novo_expr,
)

from gnomad.utils.annotations import get_copy_state_by_sex

# I want to get a table with all the following cases:
# 1. autosomal locus with HIGH P
# 2. autosomal locus with medium P
# 3. autosomal locus with low P
# 4. autosomal locus with one FAIL condition
# 5. autosomal locus with multiple FAIL conditions
# 6. hemi X locus for XY individual with a HIGH P
# 7. hemi Y locus for XY individual with a HIGH P
# 8. autosomal locus that is not de novo
# 9. autosomal locus with PLs all [0,0,0] and no freq prior
# 10. autosomal locus with missing PLs
# 11. autosomal locus with a multi-allelic site
# 12. autosomal locus with frequency prior out of range


class TestDeNovoMutation:
    """Test suite for de novo mutation functions."""

    @pytest.fixture
    def ht_de_novo_test_cases(self) -> hl.Table:
        """Fixture to create a Hail Table with different de novo mutation test cases."""
        data = [
            # 1. Autosomal locus with HIGH confidence
            {
                "locus": hl.locus("chr1", 10000, reference_genome="GRCh38"),
                "alleles": ["A", "C"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[5, 30], DP=35, GQ=99,
                                     PL=[99, 0, 1]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[20, 0], DP=20, GQ=60,
                                    PL=[0, 60, 120]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[25, 0], DP=25, GQ=80,
                                    PL=[0, 80, 150]),
                "is_xx": True,
                "freq_prior": 1e-5,
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 0.999,
                "expected_de_novo_expr": hl.struct(is_de_novo=True, p_de_novo=0.999,
                                                   confidence="HIGH",
                                                   fail_reason=hl.missing(hl.tset(hl.tstr))),
            },
            # 2. Autosomal locus with MEDIUM confidence
            {
                "locus": hl.locus("chr1", 11000, reference_genome="GRCh38"),
                "alleles": ["A", "T"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[59,61], DP=120, GQ=37,
                                     PL=[542,0,1940]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[32,0], DP=32, GQ=60,
                                    PL=[0,60,120]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[32,0], DP=37, GQ=60,
                                    PL=[0,60,120]),
                "is_xx": False,
                "freq_prior": 2.62e-03,
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 0.615,
                "expected_de_novo_expr": hl.struct(is_de_novo=True, p_de_novo=0.615,
                                                   confidence="MEDIUM",
                                                   fail_reason=hl.missing(hl.tset(hl.tstr))),
            },
            # 3. Autosomal locus with LOW confidence
            {
                "locus": hl.locus("chr1", 12000, reference_genome="GRCh38"),
                "alleles": ["G", "T"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[7,2], DP=18, GQ=43,
                                     PL=[43,0,387]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[25,0], DP=25, GQ=40,
                                    PL=[0, 40, 80]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[23,0], DP=23, GQ=40,
                                    PL=[0,40,80]),
                "is_xx": True,
                "freq_prior": 0,
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 0.926,
                "expected_de_novo_expr": hl.struct(is_de_novo=True, p_de_novo=0.926,
                                                   confidence="LOW",
                                                   fail_reason=hl.missing(
                                                       hl.tset(hl.tstr))),
            },
            # 4. Autosomal locus with one FAIL condition
            {
                "locus": hl.locus("chr1", 13000, reference_genome="GRCh38"),
                "alleles": ["C", "G"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[20, 5], DP=10, GQ=50, PL=[10, 0, 100]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[10, 0], DP=100, GQ=99, PL=[0, 99, 198]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[10, 0], DP=100, GQ=99, PL=[0, 99, 198]),
                "is_xx": True,
                "freq_prior": 1e-5,
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 1,
                "expected_de_novo_expr": hl.struct(is_de_novo=True,
                                                   p_de_novo=hl.missing(hl.tfloat64),
                                                   confidence=hl.missing(hl.tstr),
                                                   fail_reason={"min_dp_ratio"}),
            },
            # 5. Autosomal locus with multiple FAIL conditions
            {
                "locus": hl.locus("chr1", 14000, reference_genome="GRCh38"),
                "alleles": ["A", "G"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[9, 2], DP=11, GQ=2, PL=[2, 0, 230]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0, PL=[0, 0, 0]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0, PL=[0, 0, 0]),
                "is_xx": True,
                "freq_prior": 1e-5,
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 0,
                "expected_de_novo_expr": hl.struct(is_de_novo=True,
                                                   p_de_novo=hl.missing(hl.tfloat64),
                                                   confidence=hl.missing(hl.tstr),
                                                   fail_reason={"min_de_novo_p",
                            "min_proband_ab",
                            "min_proband_gq",
                            "parent_sum_ad_0"}),
            },
            # 6. Hemi X locus for XY individual with HIGH confidence
            {
                "locus": hl.locus("chrX", 8400000, reference_genome="GRCh38"),
                "alleles": ["A", "G"],
                "proband": hl.struct(GT=hl.call(1, 1), AD=[0,14], DP=14, GQ=42,
                                     PL=[419,42,0]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[38, 0], DP=38, GQ=40,
                                    PL=[0,40,80]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[97,0], DP=110, GQ=99,
                                    PL=[0,99,198]),
                "is_xx": False,
                "freq_prior": 3.74e-02,
                "expected_error": False,
                "expected_copy_state": (False, True, False),
                "expected_p_de_novo": 0.999,
                "expected_de_novo_expr": hl.struct(is_de_novo=True, p_de_novo=0.999,
                                                   confidence="HIGH",
                                                   fail_reason=hl.missing(
                                                       hl.tset(hl.tstr))),
            },
            # 7. Hemi Y locus for XY individual with HIGH confidence
            {
                "locus": hl.locus("chrY", 9900000, reference_genome="GRCh38"),
                "alleles": ["A", "G"],
                "proband": hl.struct(GT=hl.call(1, 1), AD=[0, 43], DP=43, GQ=99,
                                     PL=[1363,129,0]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[28, 0], DP=28, GQ=40,
                                    PL=[0, 40, 80]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0,
                                    PL=[0, 0, 0]),
                "is_xx": False,
                "freq_prior": hl.missing(hl.tfloat64),
                "expected_error": False,
                "expected_copy_state": (False, False, True),
                "expected_p_de_novo": 0.962,
                "expected_de_novo_expr": hl.struct(is_de_novo=True, p_de_novo=0.962,
                                                   confidence="HIGH",
                                                   fail_reason=hl.missing(
                                                       hl.tset(hl.tstr))),
            },
            # 8. Autosomal locus that is not de novo
            {
                "locus": hl.locus("chr1", 15000, reference_genome="GRCh38"),
                "alleles": ["G", "T"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[15, 10], DP=30, GQ=50, PL=[10, 0, 100]),
                "father": hl.struct(GT=hl.call(0, 1), AD=[10, 5], DP=20, GQ=40, PL=[0, 20, 80]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[20, 0], DP=20, GQ=50, PL=[0, 99, 198]),
                "is_xx": False,
                "freq_prior": 1e-5,
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 0.077,
                "expected_de_novo_expr": hl.struct(
                    is_de_novo=False,
                    p_de_novo=hl.missing(hl.tfloat64),
                    confidence=hl.missing(hl.tstr),
                    fail_reason=hl.missing(hl.tset(hl.tstr)),
                ),
            },
            # 9. Autosomal locus with PLs all [0,0,0] and no freq prior
            {
                "locus": hl.locus("chr1", 16000, reference_genome="GRCh38"),
                "alleles": ["G", "T"],
                "proband": hl.struct(GT=hl.call(0, 1), AD=[0, 2], DP=2, GQ=0,
                                     PL=[0, 0, 0]),
                "father": hl.struct(GT=hl.call(0, 0), AD=[2, 0], DP=2, GQ=0,
                                    PL=[0, 0, 0]),
                "mother": hl.struct(GT=hl.call(0, 0), AD=[2, 0], DP=2, GQ=0,
                                    PL=[0, 0, 0]),
                "is_xx": False,
                "freq_prior": hl.missing(hl.tfloat64),
                "expected_error": False,
                "expected_copy_state": (True, False, False),
                "expected_p_de_novo": 0.001,
                "expected_de_novo_expr": hl.struct(
                    is_de_novo=True,
                    p_de_novo=hl.missing(hl.tfloat64),
                    confidence=hl.missing(hl.tstr),
                    fail_reason={"min_de_novo_p",
                                 "min_proband_gq"}
                ),
            }
        ]

        # Convert list to a Hail Table
        ht = hl.Table.parallelize(
            data,
            schema=hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                proband=hl.tstruct(GT=hl.tcall, AD=hl.tarray(hl.tint32), DP=hl.tint32,
                                   GQ=hl.tint32, PL=hl.tarray(hl.tint32)),
                father=hl.tstruct(GT=hl.tcall, AD=hl.tarray(hl.tint32), DP=hl.tint32,
                                  GQ=hl.tint32, PL=hl.tarray(hl.tint32)),
                mother=hl.tstruct(GT=hl.tcall, AD=hl.tarray(hl.tint32), DP=hl.tint32,
                                  GQ=hl.tint32, PL=hl.tarray(hl.tint32)),
                is_xx=hl.tbool,
                freq_prior=hl.tfloat64,
                expected_error=hl.tbool,
                expected_copy_state=hl.ttuple(hl.tbool, hl.tbool, hl.tbool),
                expected_p_de_novo=hl.tfloat64,
                expected_de_novo_expr=hl.tstruct(
                    is_de_novo=hl.tbool,
                    p_de_novo=hl.tfloat64,
                    confidence=hl.tstr,
                    fail_reason=hl.tset(hl.tstr),
                ),
            ),
        )

        return ht

    def test_get_copy_state_by_sex(self, ht_de_novo_test_cases):
        """Test `get_copy_state_by_sex` function using a Hail Table."""
        # 🔹 Compute actual copy state using `get_copy_state_by_sex`
        ht = ht_de_novo_test_cases.annotate(
            computed_copy_state=get_copy_state_by_sex(ht_de_novo_test_cases.locus,
                                                      ht_de_novo_test_cases.is_xx)
        )

        # 🔹 Evaluate computed and expected values
        computed_values = hl.eval(ht.computed_copy_state.collect())
        expected_values = hl.eval(ht.expected_copy_state.collect())

        # 🔹 Compare expected vs. actual results
        for i, (computed, expected) in enumerate(zip(computed_values, expected_values)):
            assert computed == expected, f"Copy state mismatch at index {i}: expected {expected}, got {computed}"

    def test_calculate_de_novo_post_prob(self, ht_de_novo_test_cases):
        """Test `calculate_de_novo_post_prob` function using a Hail Table."""
        # 🔹 Store computed values and handle expected errors in Hail
        ht = ht_de_novo_test_cases.annotate(
            computed_p_de_novo=hl.case()
            .when(
                ht_de_novo_test_cases.expected_error,
                # 🔹 If error expected, return missing
                hl.missing(hl.tfloat64),
            )
            .default(
                calculate_de_novo_post_prob(
                    ht_de_novo_test_cases.proband.PL,
                    ht_de_novo_test_cases.father.PL,
                    ht_de_novo_test_cases.mother.PL,
                    ht_de_novo_test_cases.expected_copy_state[0],
                    ht_de_novo_test_cases.expected_copy_state[1],
                    ht_de_novo_test_cases.expected_copy_state[2],
                    ht_de_novo_test_cases.freq_prior,
                    min_pop_prior=100 / 3e7,
                )
            )
        )

        # 🔹 Collect the table
        ht.select("computed_p_de_novo", "expected_p_de_novo",
                            "expected_error").show(-1)
        results = ht.select("computed_p_de_novo", "expected_p_de_novo",
                            "expected_error").collect()

        for row in results:
            if row.expected_error:
                # 🔹 If an error was expected, assert the result is missing
                assert hl.is_missing(row.computed_p_de_novo)
            else:
                # 🔹 Otherwise, compare expected values
                assert round(row.computed_p_de_novo, 3) == row.expected_p_de_novo

    def test_default_get_de_novo_expr(self, ht_de_novo_test_cases):
        """Test different scenarios of `default_get_de_novo_expr` using a Hail Table."""
        # 🔹 Store computed values and handle expected errors
        ht = ht_de_novo_test_cases.annotate(
            computed_de_novo_expr=hl.case()
            .when(
                ht_de_novo_test_cases.expected_error,
                hl.missing(ht_de_novo_test_cases.expected_de_novo_expr.dtype),
            )
            .default(
                default_get_de_novo_expr(
                    ht_de_novo_test_cases.locus,
                    ht_de_novo_test_cases.alleles,
                    ht_de_novo_test_cases.proband,
                    ht_de_novo_test_cases.father,
                    ht_de_novo_test_cases.mother,
                    ht_de_novo_test_cases.is_xx,
                    ht_de_novo_test_cases.freq_prior,
                )
            )
        )

        # 🔹 Round `p_de_novo` within the struct before evaluation
        ht = ht.annotate(
            computed_de_novo_expr=hl.struct(
                is_de_novo=ht.computed_de_novo_expr.is_de_novo,
                p_de_novo=hl.or_missing(
                    hl.is_defined(ht.computed_de_novo_expr.p_de_novo),
                    hl.float64(
                        hl.int32(ht.computed_de_novo_expr.p_de_novo * 1000)) / 1000,
                ),
                confidence=ht.computed_de_novo_expr.confidence,
                fail_reason=ht.computed_de_novo_expr.fail_reason,
            )
        )

        # 🔹 Evaluate computed and expected values
        computed_values = hl.eval(ht.computed_de_novo_expr.collect())
        expected_values = hl.eval(ht.expected_de_novo_expr.collect())

        # 🔹 Compare expected vs. actual results
        for i, (computed, expected) in enumerate(zip(computed_values, expected_values)):
            assert computed == expected, f"Copy state mismatch at index {i}: expected {expected}, got {computed}"
