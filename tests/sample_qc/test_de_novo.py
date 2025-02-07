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
            ([73, 0, 161], [0, 99, 198], [0, 99, 198], True, False, False, 1.80e-02,
             100 / 3e7, 0.999),
            ([152, 0, 283], [0, 99, 198], [0, 63, 126], True, False, False, 7.55e-02,
             100 / 3e7, 0.198),
            # ❌ Invalid `freq_prior` case (should raise `HailUserError`)
            (
            [99, 50, 0], [0, 99, 198], [0, 99, 198], False, True, False, 1.2, 100 / 3e7,
            None),
        ],
    )
    def test_calculate_de_novo_post_prob(
            self, proband_pl, father_pl, mother_pl, diploid, hemi_x, hemi_y, freq_prior,
            min_pop_prior, expected
    ):
        """Test `calculate_de_novo_post_prob` function."""
        # Case where we expect an error (freq_prior is out of range)
        if expected is None:
            with pytest.raises(hl.utils.HailUserError,
                               match=r"de_novo: expect 0 <= freq_prior_expr <= 1, found .*"):
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

    def test_default_get_de_novo_expr_fail_conditions(self):
        """Test default_get_de_novo_expr with a failing case where multiple fail conditions apply."""
        # Define locus and alleles (Autosomal)
        locus = hl.locus("1", 10000)
        alleles = hl.literal(["A", "C"])

        # Define proband, father, and mother genotype structures
        proband_expr = hl.struct(
            GT=hl.call(0, 1), AD=[9, 2], DP=11, GQ=2, PL=[2, 0, 230]
        )
        father_expr = hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0, PL=[0, 0, 0])
        mother_expr = hl.struct(GT=hl.call(0, 0), AD=[0, 0], DP=0, GQ=0, PL=[0, 0, 0])

        # Population frequency prior
        freq_prior_expr = hl.literal(1e-5)
        is_xx_expr = hl.literal(True)

        # Compute de novo classification
        result_expr = default_get_de_novo_expr(
            locus,
            alleles,
            proband_expr,
            father_expr,
            mother_expr,
            is_xx_expr,
            freq_prior_expr,
        )

        # Expected result structure
        expected_result_expr = hl.struct(
            is_de_novo=True,
            p_de_novo=hl.missing(hl.tfloat64),
            confidence=hl.missing(hl.tstr),
            fail_reason=hl.set(
                {"min_de_novo_p", "min_proband_ab", "min_proband_gq", "parent_sum_ad_0"}
            ),
        )

        # Evaluate Hail expressions to convert to Python-native objects
        result = hl.eval(result_expr)
        expected_result = hl.eval(expected_result_expr)

        # Convert fail_reason to set for direct comparison
        assert result == expected_result
