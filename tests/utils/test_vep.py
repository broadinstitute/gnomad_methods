"""Tests for the VEP utility module."""

import hail as hl
import pytest

from gnomad.utils.vep import (
    get_loftee_end_trunc_filter_expr,
    update_loftee_end_trunc_filter,
)


class TestGetLofteeEndTruncFilterExpr:
    """Test the get_loftee_end_trunc_filter_expr function."""

    @pytest.fixture
    def sample_csq_structs(self):
        """Fixture to create sample consequence structs with different LOFTEE annotations."""
        return [
            # Case 1: GERP_DIST < 0, 50_BP_RULE != PASS -> should be True for default.
            hl.Struct(lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL,OTHER:value"),
            # Case 2: GERP_DIST >= 0, 50_BP_RULE != PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:1.5,50_BP_RULE:FAIL,OTHER:value"),
            # Case 3: GERP_DIST < 0, 50_BP_RULE = PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:-1.0,50_BP_RULE:PASS,OTHER:value"),
            # Case 4: GERP_DIST >= 0, 50_BP_RULE = PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:0.5,50_BP_RULE:PASS,OTHER:value"),
            # Case 5: GERP_DIST >= 0, 50_BP_RULE != PASS -> should be False for default.
            hl.Struct(lof_info="GERP_DIST:0.5,50_BP_RULE:FAIL,OTHER:value"),
            # Case 6: Missing GERP_DIST (defaults to 0), 50_BP_RULE != PASS -> should
            # be False for default.
            hl.Struct(lof_info="50_BP_RULE:FAIL,OTHER:value"),
            # Case 7: GERP_DIST < 0, missing 50_BP_RULE (defaults to empty) -> should
            # be True for default.
            hl.Struct(lof_info="GERP_DIST:-1.5,OTHER:value"),
            # Case 8: Empty lof_info -> should be False for default.
            hl.Struct(lof_info=""),
        ]

    def test_default_cutoff(self, sample_csq_structs):
        """Test the function with default cutoff of 0.0."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs],
            hl.tstruct(csq=hl.tstruct(lof_info=hl.tstr)),
        )

        # Apply the function.
        ht = ht.annotate(end_trunc=get_loftee_end_trunc_filter_expr(ht.csq))

        # Collect results
        results = ht.collect()

        # Expected results for default cutoff (0.0).
        expected = [True, False, False, False, False, False, True, False]

        assert [r.end_trunc for r in results] == expected

    def test_custom_cutoff_positive(self, sample_csq_structs):
        """Test the function with a positive cutoff."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs],
            hl.tstruct(csq=hl.tstruct(lof_info=hl.tstr)),
        )

        # Apply the function with cutoff 1.0.
        ht = ht.annotate(
            end_trunc=get_loftee_end_trunc_filter_expr(ht.csq, gerp_dist_cutoff=1.0)
        )

        # Collect results
        results = ht.collect()

        # Expected results for cutoff 1.0.
        expected = [True, False, False, False, True, True, True, True]

        assert [r.end_trunc for r in results] == expected

    def test_custom_cutoff_negative(self, sample_csq_structs):
        """Test the function with a negative cutoff."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs],
            hl.tstruct(csq=hl.tstruct(lof_info=hl.tstr)),
        )

        # Apply the function with cutoff -1.0.
        ht = ht.annotate(
            end_trunc=get_loftee_end_trunc_filter_expr(ht.csq, gerp_dist_cutoff=-1.0)
        )

        # Collect results
        results = ht.collect()

        # Expected results for cutoff -1.0.
        expected = [True, False, False, False, False, False, True, False]

        assert [r.end_trunc for r in results] == expected


class TestUpdateLofteeEndTruncFilter:
    """Test the update_loftee_end_trunc_filter function."""

    @pytest.fixture
    def sample_csq_structs_with_filters(self):
        """Fixture to create sample consequence structs with lof_filter and lof annotations."""
        return [
            # Case 1: Should add END_TRUNC filter.
            hl.Struct(
                lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL",
                lof_filter="SINGLE_EXON",
                lof="HC",
            ),
            # Case 2: Should not add END_TRUNC filter.
            hl.Struct(
                lof_info="GERP_DIST:1.5,50_BP_RULE:PASS",
                lof_filter="SINGLE_EXON",
                lof="HC",
            ),
            # Case 3: Should remove existing END_TRUNC filter.
            hl.Struct(
                lof_info="GERP_DIST:1.0,50_BP_RULE:PASS",
                lof_filter="SINGLE_EXON,END_TRUNC",
                lof="LC",
            ),
            # Case 4: Should add END_TRUNC.
            hl.Struct(
                lof_info="GERP_DIST:-1.0,50_BP_RULE:FAIL", lof_filter="", lof="HC"
            ),
            # Case 5: Missing lof_filter.
            hl.Struct(
                lof_info="GERP_DIST:-1.5,50_BP_RULE:FAIL", lof_filter=None, lof="HC"
            ),
        ]

    def test_update_single_struct(self, sample_csq_structs_with_filters):
        """Test updating a single consequence struct."""
        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in sample_csq_structs_with_filters],
            hl.tstruct(
                csq=hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csq=update_loftee_end_trunc_filter(ht.csq))

        # Collect results.
        results = ht.collect()

        # Check results.
        assert results[0].updated_csq.lof_filter == "END_TRUNC,SINGLE_EXON"
        assert results[0].updated_csq.lof == "LC"

        # Still LC because filter is not empty.
        assert results[1].updated_csq.lof_filter == "SINGLE_EXON"
        assert results[1].updated_csq.lof == "LC"

        # Still LC because filter is not empty.
        assert results[2].updated_csq.lof_filter == "SINGLE_EXON"
        assert results[2].updated_csq.lof == "LC"

        assert results[3].updated_csq.lof_filter == "END_TRUNC"
        assert results[3].updated_csq.lof == "LC"

        assert results[4].updated_csq.lof_filter == "END_TRUNC"
        assert results[4].updated_csq.lof == "LC"

    def test_update_array_of_structs(self, sample_csq_structs_with_filters):
        """Test updating an array of consequence structs."""
        # Create a table with arrays of consequences.
        ht = hl.Table.parallelize(
            [
                {"csqs": sample_csq_structs_with_filters[:2]},
                {"csqs": sample_csq_structs_with_filters[2:]},
            ],
            hl.tstruct(
                csqs=hl.tarray(
                    hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
                )
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csqs=update_loftee_end_trunc_filter(ht.csqs))

        # Collect results.
        results = ht.collect()

        # Check first array.
        first_array = results[0].updated_csqs
        assert first_array[0].lof_filter == "END_TRUNC,SINGLE_EXON"
        assert first_array[0].lof == "LC"
        assert first_array[1].lof_filter == "SINGLE_EXON"
        assert first_array[1].lof == "LC"

        # Check second array.
        second_array = results[1].updated_csqs
        assert second_array[0].lof_filter == "SINGLE_EXON"
        assert second_array[0].lof == "LC"
        assert second_array[1].lof_filter == "END_TRUNC"
        assert second_array[1].lof == "LC"
        assert second_array[2].lof_filter == "END_TRUNC"
        assert second_array[2].lof == "LC"

    def test_missing_lof_annotation(self):
        """Test updating when lof annotation is missing."""
        csq_with_missing_lof = hl.Struct(
            lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL",
            lof_filter="SINGLE_EXON",
            lof=None,
        )

        ht = hl.Table.parallelize(
            [{"csq": csq_with_missing_lof}],
            hl.tstruct(
                csq=hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csq=update_loftee_end_trunc_filter(ht.csq))

        # Collect results.
        results = ht.collect()

        # This case shouldn't happen. If lof_filter is defined, lof should be defined
        # too. However, we should handle it gracefully by adding END_TRUNC, but
        # maintaining the lof missingness status.
        assert results[0].updated_csq.lof_filter == "END_TRUNC,SINGLE_EXON"
        assert results[0].updated_csq.lof is None

    def test_empty_filter_handling(self):
        """Test handling of empty and None filters."""
        test_cases = [
            hl.Struct(
                lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL", lof_filter="", lof="HC"
            ),
            hl.Struct(
                lof_info="GERP_DIST:-2.5,50_BP_RULE:FAIL", lof_filter=None, lof="HC"
            ),
        ]

        ht = hl.Table.parallelize(
            [{"csq": csq} for csq in test_cases],
            hl.tstruct(
                csq=hl.tstruct(lof_info=hl.tstr, lof_filter=hl.tstr, lof=hl.tstr)
            ),
        )

        # Apply the function.
        ht = ht.annotate(updated_csq=update_loftee_end_trunc_filter(ht.csq))

        # Collect results.
        results = ht.collect()

        assert results[0].updated_csq.lof_filter == "END_TRUNC"
        assert results[0].updated_csq.lof == "LC"

        assert results[1].updated_csq.lof_filter == "END_TRUNC"
        assert results[1].updated_csq.lof == "LC"
