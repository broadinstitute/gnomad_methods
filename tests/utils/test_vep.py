"""Tests for the gnomAD VEP utility functions."""

from typing import Any, Dict, List, Optional

import hail as hl
import pytest

from gnomad.utils.vep import (
    add_most_severe_consequence_to_consequence,
    add_most_severe_csq_to_tc_within_vep_root,
    filter_to_most_severe_consequences,
    filter_vep_transcript_csqs_expr,
    get_most_severe_consequence_expr,
    get_most_severe_csq_from_multiple_csq_lists,
)


class TestGetMostSevereConsequenceExpr:
    """Test the get_most_severe_consequence_expr function."""

    @pytest.mark.parametrize(
        "csq_expr, custom_csq_order, expected",
        [
            # Test with default csq_order.
            (
                hl.literal(["splice_region_variant", "intron_variant"]),
                None,
                "splice_region_variant",
            ),
            # Test with custom csq_order.
            (
                hl.literal(["splice_region_variant", "intron_variant"]),
                ["intron_variant", "splice_region_variant"],
                "intron_variant",
            ),
            # Test with csq_expr that contains a consequence not in csq_order.
            (hl.literal(["non_existent_consequence"]), None, None),
            # Test with empty csq_expr.
            (hl.empty_array(hl.tstr), None, None),
        ],
    )
    def test_get_most_severe_consequence_expr(
        self,
        csq_expr: hl.expr.ArrayExpression,
        custom_csq_order: List[str],
        expected: str,
    ) -> None:
        """
        Test get_most_severe_consequence_expr with various parameters.

        :param csq_expr: The consequence terms to evaluate.
        :param custom_csq_order: The custom consequence order to use.
        :param expected: The expected most severe consequence.
        :return: None.
        """
        result = get_most_severe_consequence_expr(csq_expr, custom_csq_order)
        assert hl.eval(result) == expected, f"Expected {expected}"


class TestAddMostSevereConsequenceToConsequence:
    """Tests for the add_most_severe_consequence_to_consequence function."""

    @pytest.mark.parametrize(
        "tc_expr, custom_csq_order, expected",
        [
            # Test with default csq_order.
            (
                hl.struct(consequence_terms=["missense_variant", "synonymous_variant"]),
                None,
                "missense_variant",
            ),
            # Test with custom csq_order.
            (
                hl.struct(consequence_terms=["missense_variant", "synonymous_variant"]),
                ["synonymous_variant", "missense_variant"],
                "synonymous_variant",
            ),
            # Test with csq_expr that contains a consequence not in csq_order.
            (hl.struct(consequence_terms=["non_existent_consequence"]), None, None),
        ],
    )
    def test_add_most_severe_consequence_to_consequence(
        self,
        tc_expr: hl.expr.StructExpression,
        custom_csq_order: List[str],
        expected: str,
    ) -> None:
        """
        Test add_most_severe_consequence_to_consequence with various parameters.

        :param tc_expr: The transcript consequence to evaluate.
        :param custom_csq_order: The custom consequence order to use.
        :param expected: The expected most severe consequence.
        :return: None.
        """
        result = add_most_severe_consequence_to_consequence(tc_expr, custom_csq_order)
        assert hl.eval(result.most_severe_consequence) == hl.eval(
            get_most_severe_consequence_expr(
                tc_expr.consequence_terms, custom_csq_order
            )
        ), f"Expected {expected}"


class TestAddMostSevereCsqToTcWithinVepRoot:
    """Tests for the add_most_severe_csq_to_tc_within_vep_root function."""

    @pytest.mark.parametrize(
        "transcript_consequences, expected_most_severe",
        [
            # Test with multiple consequence terms of different severity.
            (
                [
                    hl.struct(
                        consequence_terms=["missense_variant", "synonymous_variant"]
                    )
                ],
                "missense_variant",
            ),
            # Test with a single consequence term.
            (
                [hl.struct(consequence_terms=["synonymous_variant"])],
                "synonymous_variant",
            ),
            # Test with multiple consequence terms of the same severity.
            (
                [
                    hl.struct(
                        consequence_terms=["synonymous_variant", "synonymous_variant"]
                    )
                ],
                "synonymous_variant",
            ),
            # Test with a consequence term not in the default order.
            ([hl.struct(consequence_terms=["non_existent_consequence"])], None),
        ],
    )
    def test_add_most_severe_csq_to_tc_within_vep_root(
        self,
        transcript_consequences: List[hl.expr.StructExpression],
        expected_most_severe: str,
    ) -> None:
        """
        Test the add_most_severe_csq_to_tc_within_vep_root function.

        :param transcript_consequences: List of transcript consequences to annotate.
        :param expected_most_severe: The expected most severe consequence.
        :return: None.
        """
        # Create a mock MatrixTable with vep.transcript_consequences.
        mt = hl.utils.range_matrix_table(1, 1)
        mt = mt.annotate_rows(
            vep=hl.struct(transcript_consequences=transcript_consequences)
        )

        # Apply the function.
        result_mt = add_most_severe_csq_to_tc_within_vep_root(mt)

        # Check that the most_severe_consequence is correct.
        assert (
            result_mt.vep.transcript_consequences[0].most_severe_consequence.take(1)[0]
            == expected_most_severe
        ), f"Expected '{expected_most_severe}'"


class TestFilterToMostSevereConsequences:
    """Tests for the filter_to_most_severe_consequences function."""

    @pytest.fixture
    def mock_csq_expr(self) -> hl.expr.ArrayExpression:
        """Fixture to create a mock array of VEP consequences."""

        def _build_csq_struct(
            csq: str,
            protein_coding: bool,
            lof: str,
            no_lof_flags: bool,
            polyphen_prediction: Optional[str] = None,
        ) -> hl.Struct:
            """
            Build a mock VEP consequence struct.

            :param csq: The consequence term.
            :param protein_coding: Whether the consequence is protein coding.
            :param lof: The LOF value.
            :param no_lof_flags: Whether the consequence has no LOF flags.
            :param polyphen_prediction: The PolyPhen prediction.
            :return: The mock VEP consequence struct.
            """
            return hl.Struct(
                biotype="protein_coding" if protein_coding else "not_protein_coding",
                lof=lof,
                lof_flags=hl.missing(hl.tstr) if no_lof_flags else "flag1",
                consequence_terms=[csq],
                polyphen_prediction=(
                    hl.missing(hl.tstr)
                    if polyphen_prediction is None
                    else polyphen_prediction
                ),
            )

        struct_order = [
            ["stop_gained", True, "HC", True, None],
            ["stop_lost", True, "HC", True, None],
            ["splice_acceptor_variant", True, "HC", False, "benign"],
            ["splice_acceptor_variant", True, "HC", False, "possibly_damaging"],
            ["splice_acceptor_variant", True, "LC", True, None],
            ["splice_acceptor_variant", True, "LC", True, "possibly_damaging"],
            ["splice_acceptor_variant", True, "LC", False, "possibly_damaging"],
            ["stop_gained", False, "HC", True, None],
            ["stop_gained", False, "HC", True, "probably_damaging"],
            ["splice_acceptor_variant", False, "HC", False, "possibly_damaging"],
            ["splice_acceptor_variant", False, "LC", True, "probably_damaging"],
            ["splice_acceptor_variant", False, "LC", True, None],
            ["splice_acceptor_variant", False, "LC", True, "possibly_damaging"],
            ["splice_acceptor_variant", False, "LC", False, "probably_damaging"],
        ]

        return hl.literal([_build_csq_struct(*p) for p in struct_order])

    polyphen_order = ["probably_damaging", "possibly_damaging", "benign"]
    polyphen_params = ["polyphen_prediction", polyphen_order]

    @pytest.mark.parametrize(
        "prioritize_protein_coding, prioritize_loftee, prioritize_loftee_no_flags, additional_order_field, additional_order, expected_most_severe_csq, expected_polyphen_prediction",
        [
            (False, False, False, None, None, None, None),
            (True, False, False, None, None, None, None),
            (False, True, False, None, None, None, None),
            (False, False, True, None, None, None, None),
            (False, False, False, *polyphen_params, None, None),
            (True, True, False, None, None, None, None),
            (True, False, True, None, None, None, None),
            (True, False, False, *polyphen_params, None, "possibly_damaging"),
            (False, True, True, None, None, "stop_gained", None),
            (False, True, False, *polyphen_params, None, "possibly_damaging"),
            (False, False, True, *polyphen_params, None, None),
            (True, True, True, None, None, "stop_gained", None),
            (True, True, False, *polyphen_params, None, "possibly_damaging"),
            (True, False, True, *polyphen_params, None, "possibly_damaging"),
            (False, True, True, *polyphen_params, "stop_gained", None),
            # Need to figure out class too large error
            (True, True, True, *polyphen_params, "stop_gained", "possibly_damaging"),
        ],
    )
    def test_filter_to_most_severe_consequences(
        self,
        mock_csq_expr: hl.expr.ArrayExpression,
        prioritize_protein_coding: bool,
        prioritize_loftee: bool,
        prioritize_loftee_no_flags: bool,
        additional_order_field: str,
        additional_order: List[str],
        expected_most_severe_csq: str,
        expected_polyphen_prediction: str,
    ) -> None:
        """
        Test the filter_to_most_severe_consequences function.

        :param mock_csq_expr: The mock VEP consequence expression.
        :param prioritize_protein_coding: Whether to prioritize protein coding.
        :param prioritize_loftee: Whether to prioritize LOFTEE.
        :param prioritize_loftee_no_flags: Whether to prioritize LOFTEE no flags.
        :param additional_order_field: The additional order field to use.
        :param additional_order: The additional order to use.
        :param expected_most_severe_csq: The expected most severe consequence.
        :param expected_polyphen_prediction: The expected PolyPhen prediction.
        :return: None.
        """
        result = filter_to_most_severe_consequences(
            mock_csq_expr,
            prioritize_protein_coding=prioritize_protein_coding,
            prioritize_loftee=prioritize_loftee,
            prioritize_loftee_no_flags=prioritize_loftee_no_flags,
            additional_order_field=additional_order_field,
            additional_order=additional_order,
        )

        expected_dict = hl.Struct(
            protein_coding=True,
            lof="HC",
            most_severe_consequence="splice_acceptor_variant",
            no_lof_flags=True,
        )

        def _get_csq_structs(
            csq: str,
            protein_coding: Optional[bool] = None,
            lof: Optional[str] = None,
            no_lof_flags: Optional[bool] = None,
            polyphen_prediction: Optional[str] = None,
        ) -> List[Dict[str, Any]]:
            """
            Get the expected consequence structs.

            :param csq: The consequence term to filter by.
            :param protein_coding: Whether to filter by protein coding.
            :param lof: The LOF value to filter by.
            :param no_lof_flags: Whether to filter by no LOF flags.
            :param polyphen_prediction: The PolyPhen prediction to filter by.
            :return: The expected consequence structs.
            """

            def _get_csq_criteria(s):
                keep = s.consequence_terms.contains(csq)
                if protein_coding:
                    keep &= s.biotype == "protein_coding"
                if lof:
                    keep &= s.lof == lof
                if no_lof_flags:
                    keep &= hl.is_missing(s.lof_flags)
                if polyphen_prediction:
                    keep &= s.polyphen_prediction == polyphen_prediction

                return keep

            return hl.eval(mock_csq_expr.filter(_get_csq_criteria))

        expected_polyphen_prediction = (
            expected_polyphen_prediction or "probably_damaging"
        )
        expected_most_severe_csq = expected_most_severe_csq or "splice_acceptor_variant"
        add_ms = expected_most_severe_csq != "splice_acceptor_variant"

        expected_select = (
            (["protein_coding"] if prioritize_protein_coding else [])
            + (["lof"] if prioritize_loftee else [])
            + (
                ["no_lof_flags"]
                if prioritize_loftee_no_flags or prioritize_loftee
                else []
            )
            + (["most_severe_consequence"] if not add_ms else [])
        )
        expected = expected_dict.select(
            *expected_select,
            **({"most_severe_consequence": expected_most_severe_csq} if add_ms else {}),
            consequences=_get_csq_structs(
                expected_most_severe_csq,
                protein_coding=prioritize_protein_coding,
                lof="HC" if prioritize_loftee else None,
                no_lof_flags=prioritize_loftee_no_flags,
                polyphen_prediction=(
                    expected_polyphen_prediction
                    if additional_order_field == "polyphen_prediction"
                    else None
                ),
            ),
        )
        assert hl.eval(result) == expected, f"Expected '{expected}'"


class TestFilterVepTranscriptCsqsExprLoftee:
    """Tests for the filter_vep_transcript_csqs_expr function."""

    @pytest.fixture
    def csq_expr(self) -> hl.expr.ArrayExpression:
        """Fixture to create a mock array of VEP consequences with LOFTEE annotations."""
        return hl.literal(
            [
                hl.struct(lof="HC", lof_flags=hl.missing(hl.tstr)),
                hl.struct(lof="HC", lof_flags=""),
                hl.struct(lof="LC", lof_flags="flag1"),
                hl.struct(lof="OS", lof_flags=hl.missing(hl.tstr)),
                hl.struct(lof=hl.missing(hl.tstr), lof_flags="flag2"),
            ]
        )

    @staticmethod
    def check_length(result: hl.expr.ArrayExpression, expected_len: int) -> None:
        """
        Check the length of the result.

        :param result: The result to check.
        :param expected_len: The expected length.
        :return: None.
        """
        assert (
            hl.eval(hl.len(result)) == expected_len
        ), f"Expected {expected_len} consequences"

    @staticmethod
    def check_lof(result: hl.expr.ArrayExpression, expected_lof: str) -> None:
        """
        Check the LOF annotation.

        :param result: The result to check.
        :param expected_lof: The expected LOF annotation.
        :return: None.
        """
        assert hl.eval(result[0].lof) == expected_lof, f"Expected '{expected_lof}'"

    @staticmethod
    def check_no_lof_flags(result: hl.expr.ArrayExpression) -> None:
        """
        Check the no LOF flags value.

        :param result: The result to check.
        :return: None.
        """
        assert hl.eval(
            hl.all(
                result.map(
                    lambda csq: hl.is_missing(csq.lof_flags) | (csq.lof_flags == "")
                )
            )
        ), "Expected no LOFTEE flags"

    @pytest.mark.parametrize(
        "loftee_labels, no_lof_flags, expected_len, expected_lof, expected_no_lof_flags",
        [
            (None, None, 5, None, None),
            (["HC"], None, 2, "HC", None),
            (None, True, 3, None, True),
            (["HC"], True, 2, "HC", True),
        ],
    )
    def test_filter_vep_transcript_csqs_expr_loftee(
        self,
        csq_expr: hl.expr.ArrayExpression,
        loftee_labels: list,
        no_lof_flags: bool,
        expected_len: int,
        expected_lof: str,
        expected_no_lof_flags: bool,
    ) -> None:
        """
        Test the filter_vep_transcript_csqs_expr function.

        :param csq_expr: The VEP consequence expression to filter.
        :param loftee_labels: The LOFTEE labels to filter by.
        :param no_lof_flags: Whether to filter by no LOF flags.
        :param expected_len: The expected length of the result.
        :param expected_lof: The expected LOF value.
        :param expected_no_lof_flags: The expected no LOF flags value.
        :return: None.
        """
        result = filter_vep_transcript_csqs_expr(
            csq_expr, loftee_labels=loftee_labels, no_lof_flags=no_lof_flags
        )
        self.check_length(result, expected_len)
        if expected_lof is not None:
            self.check_lof(result, expected_lof)
        if expected_no_lof_flags is not None:
            self.check_no_lof_flags(result)


class TestGetMostSevereCsqFromMultipleCsqLists:
    """Tests for the get_most_severe_csq_from_multiple_csq_lists function."""

    @pytest.fixture
    def vep_expr(self) -> hl.expr.StructExpression:
        """Fixture to create a mock VEP expression."""
        return hl.struct(
            transcript_consequences=[
                hl.struct(
                    biotype="protein_coding",
                    lof="HC",
                    lof_flags="flag1",
                    consequence_terms=["splice_acceptor_variant"],
                ),
                hl.struct(
                    biotype="protein_coding",
                    lof="HC",
                    lof_flags="flag1",
                    consequence_terms=["splice_acceptor_variant"],
                ),
                hl.struct(
                    biotype="protein_coding",
                    lof="HC",
                    lof_flags=hl.missing(hl.tstr),
                    consequence_terms=["stop_lost"],
                ),
                hl.struct(
                    biotype="protein_coding",
                    lof="HC",
                    lof_flags=hl.missing(hl.tstr),
                    consequence_terms=["stop_gained"],
                ),
                hl.struct(
                    biotype="protein_coding",
                    lof=hl.missing(hl.tstr),
                    lof_flags=hl.missing(hl.tstr),
                    consequence_terms=["missense_variant"],
                ),
            ],
            intergenic_consequences=[
                hl.struct(
                    biotype="intergenic",
                    consequence_terms=["intergenic_variant"],
                )
            ],
        )

    @pytest.mark.parametrize(
        "prioritize_loftee_no_flags, include_csqs, expected_most_severe, expected_protein_coding, expected_lof, expected_no_lof_flags, expected_transcript_consequences_len",
        [
            (False, False, "splice_acceptor_variant", True, "HC", True, None),
            (True, False, "stop_gained", True, "HC", True, None),
            (False, True, "splice_acceptor_variant", True, "HC", True, 2),
            (True, True, "stop_gained", True, "HC", True, 1),
        ],
    )
    def test_get_most_severe_csq_from_multiple_csq_lists(
        self,
        vep_expr: hl.expr.StructExpression,
        prioritize_loftee_no_flags: bool,
        include_csqs: bool,
        expected_most_severe: str,
        expected_protein_coding: bool,
        expected_lof: str,
        expected_no_lof_flags: bool,
        expected_transcript_consequences_len: Optional[int],
    ) -> None:
        """
        Test the get_most_severe_csq_from_multiple_csq_lists function.

        :param vep_expr: The VEP expression to test.
        :param prioritize_loftee_no_flags: Whether to prioritize LOFTEE no flags.
        :param include_csqs: Whether to include consequences.
        :param expected_most_severe: The expected most severe consequence.
        :param expected_protein_coding: The expected protein coding value.
        :param expected_lof: The expected LOF value.
        :param expected_no_lof_flags: The expected no LOF flags value.
        :param expected_transcript_consequences_len: The expected length of transcript
            consequences.
        :return: None.
        """
        result = get_most_severe_csq_from_multiple_csq_lists(
            vep_expr,
            prioritize_loftee_no_flags=prioritize_loftee_no_flags,
            include_csqs=include_csqs,
        )
        self.check_most_severe_consequence(result, expected_most_severe)
        self.check_protein_coding(result, expected_protein_coding)
        self.check_lof(result, expected_lof)
        self.check_no_lof_flags(result, expected_no_lof_flags)
        if include_csqs:
            self.check_transcript_consequences_len(
                result, expected_transcript_consequences_len
            )

    @staticmethod
    def check_most_severe_consequence(
        result: hl.expr.StructExpression, expected: str
    ) -> None:
        """
        Check the most severe consequence.

        :param result: The result to check.
        :param expected: The expected most severe consequence.
        :return: None.
        """
        assert (
            hl.eval(result.most_severe_consequence) == expected
        ), f"Expected '{expected}'"

    @staticmethod
    def check_protein_coding(result: hl.expr.StructExpression, expected: bool) -> None:
        """
        Check the protein coding value.

        :param result: The result to check.
        :param expected: The expected protein coding value.
        :return: None.
        """
        assert hl.eval(result.protein_coding) == expected, f"Expected '{expected}'"

    @staticmethod
    def check_lof(result: hl.expr.StructExpression, expected: str) -> None:
        """
        Check the LOF value.

        :param result: The result to check.
        :param expected: The expected LOF value.
        :return: None.
        """
        assert hl.eval(result.lof) == expected, f"Expected '{expected}'"

    @staticmethod
    def check_no_lof_flags(result: hl.expr.StructExpression, expected: bool) -> None:
        """
        Check the no LOF flags value.

        :param result: The result to check.
        :param expected: The expected no LOF flags value.
        :return: None.
        """
        assert hl.eval(result.no_lof_flags) == expected, f"Expected '{expected}'"

    @staticmethod
    def check_transcript_consequences_len(
        result: hl.expr.StructExpression, expected: int
    ) -> None:
        """
        Check the length of transcript consequences.

        :param result: The result to check.
        :param expected: The expected length.
        :return: None.
        """
        assert (
            hl.eval(hl.len(result.transcript_consequences)) == expected
        ), f"Expected {expected} transcript consequences"
