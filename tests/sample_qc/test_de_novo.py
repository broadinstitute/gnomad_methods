"""Test suite for de novo mutation functions."""

import hail as hl
import pytest

from gnomad.sample_qc.relatedness import (
    calculate_de_novo_post_prob,
    default_get_de_novo_expr,
)


class TestDeNovoMutation:
    """Test suite for de novo mutation functions."""

    @pytest.mark.parametrize(
        "proband_pl, father_pl, mother_pl, diploid, hemi_x, hemi_y, freq_prior, min_pop_prior, expected",
        [
            # ✅ Valid test cases (should return numeric values)
            (
                [73, 0, 161],
                [0, 99, 198],
                [0, 99, 198],
                True,
                False,
                False,
                1.80e-02,
                100 / 3e7,
                0.999,
            ),
            (
                [152, 0, 283],
                [0, 99, 198],
                [0, 63, 126],
                True,
                False,
                False,
                7.55e-02,
                100 / 3e7,
                0.198,
            ),
            # ❌ Invalid `freq_prior` case (should raise `HailUserError`)
            (
                [99, 50, 0],
                [0, 99, 198],
                [0, 99, 198],
                False,
                True,
                False,
                1.2,
                100 / 3e7,
                None,
            ),
        ],
    )
    def test_calculate_de_novo_post_prob(
        self,
        proband_pl,
        father_pl,
        mother_pl,
        diploid,
        hemi_x,
        hemi_y,
        freq_prior,
        min_pop_prior,
        expected,
    ):
        """Test `calculate_de_novo_post_prob` function."""
        # Case where we expect an error (freq_prior is out of range)
        if expected is None:
            with pytest.raises(
                hl.utils.HailUserError,
                match=r"de_novo: expect 0 <= freq_prior_expr <= 1, found .*",
            ):
                hl.eval(
                    calculate_de_novo_post_prob(
                        hl.literal(proband_pl),
                        hl.literal(father_pl),
                        hl.literal(mother_pl),
                        hl.literal(diploid),
                        hl.literal(hemi_x),
                        hl.literal(hemi_y),
                        hl.literal(freq_prior),  # Invalid frequency prior
                        min_pop_prior,
                    )
                )
        else:
            # Case where we expect a valid float result
            p_dn_expr = calculate_de_novo_post_prob(
                hl.literal(proband_pl),
                hl.literal(father_pl),
                hl.literal(mother_pl),
                hl.literal(diploid),
                hl.literal(hemi_x),
                hl.literal(hemi_y),
                hl.literal(freq_prior),
                min_pop_prior,
            )

            # Assert with floating-point tolerance
            assert round(hl.eval(p_dn_expr), 3) == expected

    @pytest.mark.parametrize(
        "locus, alleles, proband_expr, father_expr, mother_expr, is_xx_expr, freq_prior_expr, expected",
        [
            # Case 1: Multiple fail conditions
            (
                hl.locus("chr1", 10000, reference_genome="GRCh38"),
                hl.literal(["A", "C"]),
                hl.struct(GT=hl.call(0, 1), AD=[9, 2], DP=11, GQ=2, PL=[2, 0, 230]),
                hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0, PL=[0, 0, 0]),
                hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0, PL=[0, 0, 0]),
                hl.literal(True),
                hl.literal(1e-5),
                hl.struct(
                    is_de_novo=True,
                    p_de_novo=hl.missing(hl.tfloat64),
                    confidence=hl.missing(hl.tstr),
                    fail_reason=hl.set(
                        {
                            "min_de_novo_p",
                            "min_proband_ab",
                            "min_proband_gq",
                            "parent_sum_ad_0",
                        }
                    ),
                ),
            ),
            # Case 2: One fail condition (low DP ratio)
            (
                hl.locus("chr1", 20000, reference_genome="GRCh38"),
                hl.literal(["A", "T"]),
                hl.struct(GT=hl.call(0, 1), AD=[20, 5], DP=10, GQ=50, PL=[10, 0, 100]),
                hl.struct(GT=hl.call(0, 0), AD=[10, 0], DP=100, GQ=99, PL=[0, 99, 198]),
                hl.struct(GT=hl.call(0, 0), AD=[10, 0], DP=100, GQ=99, PL=[0, 99, 198]),
                hl.literal(False),
                hl.literal(1e-5),
                hl.struct(
                    is_de_novo=True,
                    p_de_novo=hl.missing(hl.tfloat64),
                    confidence=hl.missing(hl.tstr),
                    fail_reason=hl.set({"min_dp_ratio"}),
                ),
            ),
            # Case 3: Variant is inherited (not de novo)
            (
                hl.locus("chr1", 30000, reference_genome="GRCh38"),
                hl.literal(["G", "T"]),
                hl.struct(GT=hl.call(0, 1), AD=[15, 10], DP=30, GQ=50, PL=[10, 0, 100]),
                hl.struct(GT=hl.call(0, 1), AD=[10, 5], DP=20, GQ=40, PL=[0, 20, 80]),
                hl.struct(GT=hl.call(0, 0), AD=[20, 0], DP=20, GQ=50, PL=[0, 99, 198]),
                hl.literal(True),
                hl.literal(1e-5),
                hl.struct(
                    is_de_novo=False,
                    p_de_novo=hl.missing(hl.tfloat64),
                    confidence=hl.missing(hl.tstr),
                    fail_reason=hl.empty_set(hl.tstr),
                ),
            ),
            # Case 4: A passing case (high confidence de novo)
            (
                hl.locus("chr1", 40000, reference_genome="GRCh38"),
                hl.literal(["C", "G"]),
                hl.struct(GT=hl.call(0, 1), AD=[5, 30], DP=35, GQ=99, PL=[99, 0, 1]),
                hl.struct(GT=hl.call(0, 0), AD=[20, 0], DP=20, GQ=60, PL=[0, 60, 120]),
                hl.struct(GT=hl.call(0, 0), AD=[25, 0], DP=25, GQ=80, PL=[0, 80, 150]),
                hl.literal(True),
                hl.literal(1e-5),
                hl.struct(
                    is_de_novo=True,
                    p_de_novo=0.999,  # High confidence P(de novo)
                    confidence="HIGH",
                    fail_reason=hl.empty_set(hl.tstr),
                ),
            ),
        ],
    )
    def test_default_get_de_novo_expr(
        self,
        locus,
        alleles,
        proband_expr,
        father_expr,
        mother_expr,
        is_xx_expr,
        freq_prior_expr,
        expected,
    ):
        """Test different scenarios of `default_get_de_novo_expr` in one function."""
        result_expr = default_get_de_novo_expr(
            locus,
            alleles,
            proband_expr,
            father_expr,
            mother_expr,
            is_xx_expr,
            freq_prior_expr,
        )

        result = hl.eval(result_expr)
        expected_result = hl.eval(expected)

        assert result.is_de_novo == expected_result.is_de_novo
        assert (
            None if result.p_de_novo is None else round(result.p_de_novo, 3)
        ) == expected_result.p_de_novo
        assert result.confidence == expected_result.confidence
        assert result.fail_reason == expected_result.fail_reason
