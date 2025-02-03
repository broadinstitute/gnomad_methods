"""Test suite for de novo mutation functions."""

import pytest
import hail as hl

from gnomad.sample_qc.relatedness import (
    get_freq_prior,
    transform_pl_to_pp,
    calculate_de_novo_post_prob,
    call_de_novo,
    get_de_novo_expr,
)

from gnomad.utils.annotations import get_copy_state_by_sex


class TestDeNovoMutation:
    """Test suite for de novo mutation functions."""

    loci: dict[str, hl.expr.LocusExpression]

    @classmethod
    def setup_class(cls):
        """Set up common test data for all tests."""
        cls.locus_expr = hl.locus("1", 123456)
        cls.alleles_expr = hl.literal(["A", "T"])
        cls.freq_prior_expr = hl.literal(0.01)
        cls.is_xx_expr = hl.literal(False)

        # Mock Genotype Likelihoods (PL)
        cls.proband_pl = hl.literal([0, 10, 100])
        cls.father_pl = hl.literal([0, 100, 100])
        cls.mother_pl = hl.literal([0, 100, 100])

        # Mock Genotype Calls
        cls.proband_expr = hl.struct(
            GT=hl.call(0, 1), DP=10, GQ=30, AD=[3, 7], PL=cls.proband_pl
        )
        cls.father_expr = hl.struct(
            GT=hl.call(0, 0), DP=12, GQ=40, AD=[12, 0], PL=cls.father_pl
        )
        cls.mother_expr = hl.struct(
            GT=hl.call(0, 0), DP=15, GQ=50, AD=[15, 0], PL=cls.mother_pl
        )

        cls.loci = {
            "autosomal": hl.locus("chr1", 100000, reference_genome="GRCh38"),
            # PAR regions (always diploid)
            "par1": hl.locus("chrX", 2781479, reference_genome="GRCh38"),  # PAR1 start
            "par2": hl.locus("chrX", 155701383, reference_genome="GRCh38"),  # PAR2 end
            # X non-PAR (diploid for XX, hemizygous for XY)
            "x_nonpar": hl.locus("chrX", 3000000, reference_genome="GRCh38"),
            # Y non-PAR (hemizygous for XY)
            "y_nonpar": hl.locus("chrY", 10000000, reference_genome="GRCh38"),
        }

    @pytest.mark.parametrize(
        "freq_prior, min_pop_prior, expect_error, expected",
        [
            (0.05, 100 / 3e7, False, 0.05),
            (-0.01, 100 / 3e7, True, None),
            (1.2, 100 / 3e7, True, None),
            (hl.missing(hl.tfloat64), 100 / 3e7, False, 100 / 3e7),
        ],
    )
    def test_get_freq_prior(
        self, freq_prior, min_pop_prior, expect_error, expected
    ) -> None:
        """
        Test frequency prior computation.

        :param freq_prior: Frequency prior value.
        :param min_pop_prior: Minimum population prior.
        :param expect_error: Whether an error is expected.
        :param expected: Expected frequency prior.
        :return: None.
        """
        if expect_error:
            with pytest.raises(
                hl.utils.java.HailUserError,
                match="de_novo: expect 0 <= freq_prior_expr <= 1",
            ):
                expr = get_freq_prior(hl.literal(freq_prior), min_pop_prior)
                hl.eval(expr)  # Hail will throw an error here
        else:
            expr = get_freq_prior(hl.literal(freq_prior), min_pop_prior)
            result = hl.eval(expr)  # Evaluate the expression
            assert result == pytest.approx(
                expected, rel=1e-6
            )  # Compare floating point values safely

    @pytest.mark.parametrize(
        "pl_input, expected",
        [
            (
                [0, 10, 100],
                [0.9090909090082644, 0.09090909090082644, 9.090909090082645e-11],
            ),
            ([0, 0, 0], [0.3333333333333333, 0.3333333333333333, 0.3333333333333333]),
        ],
    )
    def test_transform_pl_to_pp(self, pl_input, expected) -> None:
        """
        Test PL to PP transformation.

        :param pl_input: Input PL values.
        :param expected: Expected PP values.
        :return: None.
        """
        expr = transform_pl_to_pp(hl.literal(pl_input))
        result = hl.eval(expr)

        assert result == pytest.approx(
            expected, abs=1e-12
        ), f"Got {result}, expected {expected}"

    @pytest.mark.parametrize(
        "locus_key, is_xx, expected_diploid, expected_hemi_x, expected_hemi_y",
        [
            ("autosomal", True, True, False, False),
            ("autosomal", False, True, False, False),
            ("par1", True, True, False, False),
            ("par2", False, True, False, False),
            ("x_nonpar", True, True, False, False),
            ("x_nonpar", False, False, True, False),
            ("y_nonpar", True, False, False, False),
            ("y_nonpar", False, False, False, True),
        ],
    )
    def test_get_copy_state_by_sex(
        self, locus_key, is_xx, expected_diploid, expected_hemi_x, expected_hemi_y
    ) -> None:
        """
        Test copy state determination based on locus type and sex.

        :param locus_key: Locus key.
        :param is_xx: Whether the individual is XX.
        :param expected_diploid: Expected diploid state.
        :param expected_hemi_x: Expected hemizygous X state.
        :param expected_hemi_y: Expected hemizygous Y state.
        :return: None.
        """
        locus = self.loci[locus_key]
        is_xx_expr = hl.literal(is_xx)

        diploid, hemi_x, hemi_y = get_copy_state_by_sex(locus, is_xx_expr)
        result = hl.eval([diploid, hemi_x, hemi_y])

        assert result == [
            expected_diploid,
            expected_hemi_x,
            expected_hemi_y,
        ], f"Failed for locus={locus}, is_xx={is_xx}. Expected {[expected_diploid, expected_hemi_x, expected_hemi_y]}, got {result}"

    @pytest.mark.parametrize(
        "locus_key, proband_gt, father_gt, mother_gt, is_xx, expected",
        [
            ("autosomal", (0, 1), (0, 0), (0, 0), False, True),
            ("autosomal", (1, 1), (0, 0), (0, 0), False, False),
            ("x_nonpar", (1, 1), None, (0, 0), False, True),
            ("x_nonpar", (1, 1), (0, 0), None, False, None),
            ("y_nonpar", (1, 1), None, None, False, None),
        ],
    )
    def test_call_de_novo(
        self, locus_key, proband_gt, father_gt, mother_gt, is_xx, expected
    ) -> None:
        """
        Test de novo mutation detection with different loci and parental genotypes.

        :param locus_key: Locus key.
        :param proband_gt: Proband genotype.
        :param father_gt: Father genotype.
        :param mother_gt: Mother genotype.
        :param is_xx: Whether the individual is XX.
        :param expected: Expected de novo mutation status.
        :return: None.
        """
        locus_expr = self.loci[locus_key]
        proband_expr = hl.struct(
            GT=hl.call(*proband_gt) if proband_gt else hl.missing(hl.tcall)
        )
        father_expr = hl.struct(
            GT=hl.call(*father_gt) if father_gt else hl.missing(hl.tcall)
        )
        mother_expr = hl.struct(
            GT=hl.call(*mother_gt) if mother_gt else hl.missing(hl.tcall)
        )
        is_xx_expr = hl.literal(is_xx)

        expr = call_de_novo(
            locus_expr, proband_expr, father_expr, mother_expr, is_xx_expr
        )
        result = hl.eval(expr)

        assert (
            result == expected
        ), f"Mismatch in {locus_key}: Expected {expected}, got {result}"

    def test_calculate_de_novo_post_prob(self):
        """Test posterior probability computation for de novo mutations."""
        expr = calculate_de_novo_post_prob(
            self.proband_pl,
            self.father_pl,
            self.mother_pl,
            self.locus_expr,
            self.is_xx_expr,
            self.freq_prior_expr,
        )
        result = hl.eval(expr)
        assert 0 <= result <= 1  # Posterior probability should be within valid range

    def test_get_de_novo_expr(self):
        """Test the de novo expression struct output."""
        expr = get_de_novo_expr(
            self.locus_expr,
            self.alleles_expr,
            self.proband_expr,
            self.father_expr,
            self.mother_expr,
            self.is_xx_expr,
            self.freq_prior_expr,
        )
        result = hl.eval(expr)

        assert "p_de_novo" in result
        assert "confidence" in result
        assert 0 <= result.p_de_novo <= 1  # Probability must be valid
