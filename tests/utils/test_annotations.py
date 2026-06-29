"""Tests for the annotations utility module."""

from typing import Dict, List

import ga4gh.core as ga4gh_core
import ga4gh.vrs as ga4gh_vrs
import hail as hl
import pytest

from gnomad.utils.annotations import (
    VRS_CHROM_IDS,
    add_gks_va,
    add_gks_vrs,
    annotate_downsamplings,
    annotate_freq,
    check_annotation_missingness,
    expand_strata_array_from_leaves,
    fill_missing_key_combinations,
    find_minimal_strata_groups,
    generate_freq_group_membership_array,
    get_copy_state_by_sex,
    merge_array_expressions,
    merge_freq_arrays,
    merge_histograms,
    missing_struct_expr,
)


class TestFillMissingKeyCombinations:
    """Test the fill_missing_key_combinations function."""

    expected_result = [
        hl.Struct(key1="A", key2=1, value=10),
        hl.Struct(key1="A", key2=2, value=20),
        hl.Struct(key1="B", key2=1, value=30),
        hl.Struct(key1="B", key2=2, value=None),
    ]

    expected_result_extra = [
        hl.Struct(key1="A", key2=1, value=10),
        hl.Struct(key1="A", key2=2, value=20),
        hl.Struct(key1="B", key2=1, value=30),
        hl.Struct(key1="B", key2=2, value=None),
        hl.Struct(key1="C", key2=1, value=None),
        hl.Struct(key1="C", key2=2, value=None),
    ]

    @pytest.fixture
    def sample_hail_table(self):
        """Fixture to create a sample Hail Table."""
        return hl.Table.parallelize(
            [
                {"key1": "A", "key2": 1, "value": 10},
                {"key1": "A", "key2": 2, "value": 20},
                {"key1": "B", "key2": 1, "value": 30},
            ],
            hl.tstruct(key1=hl.tstr, key2=hl.tint32, value=hl.tint32),
        ).key_by("key1", "key2")

    @pytest.mark.parametrize(
        "expected,key_values",
        [
            (expected_result, None),
            (expected_result, {"key1": ["A", "B"], "key2": [1, 2]}),
            (expected_result_extra, {"key1": ["A", "B", "C"], "key2": [1, 2]}),
            (
                expected_result_extra,
                {"key1": ["A", "B", "C"], "key2": [1, 2], "key3": ["X"]},
            ),
        ],
    )
    def test_fill_missing_key_combinations(
        self,
        sample_hail_table: hl.Table,
        expected: List[hl.Struct],
        key_values: Dict[str, List[str]],
    ) -> None:
        """Test the `fill_missing_key_combinations` function."""
        # Define fill values.
        fill_values = {"value": hl.missing(hl.tint32)}

        # Call the function.
        result_ht = fill_missing_key_combinations(
            sample_hail_table, fill_values, key_values
        )

        # Collect results.
        result = result_ht.collect()

        # Verify the result.
        assert result == expected


def test_missing_struct_expr() -> None:
    """Test the missing_struct_expr function."""
    # Define a sample tstruct.
    dtypes = hl.tstruct(field1=hl.tint32, field2=hl.tstr, field3=hl.tfloat64)

    # Call the function.
    result = missing_struct_expr(dtypes)

    # Expected result.
    expected = hl.struct(
        field1=hl.missing(hl.tint32),
        field2=hl.missing(hl.tstr),
        field3=hl.missing(hl.tfloat64),
    )

    # Verify the result.
    assert hl.eval(result == expected)


class TestGetCopyStateBySex:
    """Test the `get_copy_state_by_sex` function."""

    @pytest.mark.parametrize(
        "locus, is_xx, expected_diploid, expected_hemi_x, expected_hemi_y",
        [
            (
                hl.locus("chr1", 100000, reference_genome="GRCh38"),
                True,
                True,
                False,
                False,
            ),
            (
                hl.locus("chrX", 2781479, reference_genome="GRCh38"),
                False,
                True,
                False,
                False,
            ),
            (
                hl.locus("chrX", 3000000, reference_genome="GRCh38"),
                False,
                False,
                True,
                False,
            ),
            (
                hl.locus("chrY", 10000000, reference_genome="GRCh38"),
                False,
                False,
                False,
                True,
            ),
        ],
    )
    def test_get_copy_state_by_sex(
        self, locus, is_xx, expected_diploid, expected_hemi_x, expected_hemi_y
    ) -> None:
        """Test copy state determination based on locus type and sex."""
        is_xx_expr = hl.literal(is_xx)

        diploid, hemi_x, hemi_y = get_copy_state_by_sex(locus, is_xx_expr)
        result = hl.eval([diploid, hemi_x, hemi_y])

        assert result == [
            expected_diploid,
            expected_hemi_x,
            expected_hemi_y,
        ], f"Failed for locus={locus}, is_xx={is_xx}. Expected {[expected_diploid, expected_hemi_x, expected_hemi_y]}, got {result}"


class TestMergeArrayExpressions:
    """Test the merge_array_expressions function."""

    @pytest.fixture
    def sample_ht(self):
        """Create a sample Hail Table for testing."""
        return hl.Table.parallelize(
            [
                {"id": 1, "group": "A", "value1": 10, "value2": 20},
                {"id": 2, "group": "B", "value1": 30, "value2": 40},
                {"id": 3, "group": "A", "value1": 50, "value2": 60},
            ],
            hl.tstruct(id=hl.tint32, group=hl.tstr, value1=hl.tint32, value2=hl.tint32),
        )

    def test_merge_integer_arrays_sum(self, sample_ht):
        """Test merging integer arrays with sum operation."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20, 30]),
            array2=hl.array([5, 15, 25]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        meta2 = [{"group": "A"}, {"group": "B"}, {"group": "D"}]

        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2], [meta1, meta2], operation="sum"
        )[:2]

        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}, {"group": "D"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 4
            # Check actual values: A: 10+5=15, B: 20+15=35, C: 30+0=30, D: 0+25=25.
            assert row.result_array[0] == 15  # Group A.
            assert row.result_array[1] == 35  # Group B.
            assert row.result_array[2] == 30  # Group C.
            assert row.result_array[3] == 25  # Group D.

    def test_merge_integer_arrays_diff(self, sample_ht):
        """Test merging integer arrays with diff operation."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20, 30]),
            array2=hl.array([5, 15, 25]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        meta2 = [{"group": "A"}, {"group": "B"}, {"group": "D"}]

        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2], [meta1, meta2], operation="diff"
        )[:2]

        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 3
            # Check actual values: A: 10-5=5, B: 20-15=5, C: 30-0=30.
            assert row.result_array[0] == 5  # Group A.
            assert row.result_array[1] == 5  # Group B.
            assert row.result_array[2] == 30  # Group C.

    def test_merge_struct_arrays_sum(self, sample_ht):
        """Test merging struct arrays with sum operation."""
        ht = sample_ht.annotate(
            array1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            array2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            struct_fields=["AC", "AN", "homozygote_count"],
        )[:2]

        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 3
            # Check first element (group A): AC=10+5=15, AN=100+50=150,
            # homozygote_count=2+1=3.
            first_struct = row.result_array[0]
            assert first_struct.AC == 15
            assert first_struct.AN == 150
            assert first_struct.homozygote_count == 3

    def test_merge_struct_arrays_diff(self, sample_ht):
        """Test merging struct arrays with diff operation."""
        ht = sample_ht.annotate(
            array1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            array2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="diff",
            struct_fields=["AC", "AN", "homozygote_count"],
        )[:2]

        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 2
            # Check first element (group A): AC=10-5=5, AN=100-50=50,
            # homozygote_count=2-1=1.
            first_struct = row.result_array[0]
            assert first_struct.AC == 5
            assert first_struct.AN == 50
            assert first_struct.homozygote_count == 1

    def test_merge_with_negatives_error(self, sample_ht):
        """Test that negative values raise error when set_negatives_to_zero=False."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([15, 25]),  # Will cause negative values in diff operation.
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Get the result array (this won't raise an error yet).
        result_array = merge_array_expressions(
            [ht.array1, ht.array2], [meta1, meta2], operation="diff"
        )[0]

        # The error should be raised when we try to evaluate the result.
        with pytest.raises(Exception):
            ht.select(result_array=result_array).collect()

    def test_merge_with_negatives_set_to_zero(self, sample_ht):
        """Test that negative values are set to zero when set_negatives_to_zero=True."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([15, 25]),  # Will cause negative values in diff operation.
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Should not raise error and set negatives to zero.
        result_array = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="diff",
            set_negatives_to_zero=True,
        )[0]

        result = ht.select(result_array=result_array).collect()

        # Check that negative values are set to zero.
        for row in result:
            # A: 10-15=-5 -> 0, B: 20-25=-5 -> 0.
            assert row.result_array[0] == 0
            assert row.result_array[1] == 0

    def test_invalid_operation(self, sample_ht):
        """Test that invalid operation raises ValueError."""
        ht = sample_ht.annotate(array1=hl.array([1, 2]), array2=hl.array([3, 4]))
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        with pytest.raises(
            ValueError, match="Operation must be either 'sum' or 'diff'!"
        ):
            merge_array_expressions(
                [ht.array1, ht.array2], [meta1, meta2], operation="invalid"
            )

    def test_insufficient_arrays(self, sample_ht):
        """Test that insufficient arrays raises ValueError."""
        ht = sample_ht.annotate(array1=hl.array([1, 2]))
        meta1 = [{"group": "A"}, {"group": "B"}]

        with pytest.raises(
            ValueError, match="Must provide at least two arrays to merge!"
        ):
            merge_array_expressions([ht.array1], [meta1])

    def test_mismatched_lengths(self, sample_ht):
        """Test that mismatched array and meta lengths raises ValueError."""
        ht = sample_ht.annotate(array1=hl.array([1, 2]), array2=hl.array([3, 4]))
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]
        meta3 = [{"group": "A"}]  # Different length

        with pytest.raises(
            ValueError, match="Length of arrays and meta must be equal!"
        ):
            merge_array_expressions([ht.array1, ht.array2], [meta1, meta2, meta3])

    def test_merge_array_expressions_with_count_arrays(self, sample_ht):
        """Test merging arrays with count_arrays argument."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
        )

        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 3
            assert len(row.result_counts["sample_count"]) == 3
            # Check main array: A: 10+5=15, B: 20+0=20, C: 0+15=15
            assert row.result_array[0] == 15
            assert row.result_array[1] == 20
            assert row.result_array[2] == 15
            # Check count array: A: 100+50=150, B: 200+0=200, C: 0+150=150
            assert row.result_counts["sample_count"][0] == 150
            assert row.result_counts["sample_count"][1] == 200
            assert row.result_counts["sample_count"][2] == 150

    def test_merge_struct_arrays_with_count_arrays(self, sample_ht):
        """Test merging struct arrays with count_arrays argument."""
        ht = sample_ht.annotate(
            array1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            array2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            struct_fields=["AC", "AN", "homozygote_count"],
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 3
            assert len(row.result_counts["sample_count"]) == 3
            # Check first struct (group A): AC=10+5=15, AN=100+50=150,
            # homozygote_count=2+1=3
            first_struct = row.result_array[0]
            assert first_struct.AC == 15
            assert first_struct.AN == 150
            assert first_struct.homozygote_count == 3
            # Check count array: A: 100+50=150, B: 200+0=200, C: 0+150=150
            assert row.result_counts["sample_count"][0] == 150
            assert row.result_counts["sample_count"][1] == 200
            assert row.result_counts["sample_count"][2] == 150

    def test_merge_arrays_with_multiple_count_arrays(self, sample_ht):
        """Test merging arrays with multiple count_arrays."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
        )

        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
            qual1=hl.array([30, 40]),
            qual2=hl.array([25, 35]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            count_arrays={
                "sample_count": [ht.count1, ht.count2],
                "quality_score": [ht.qual1, ht.qual2],
            },
        )

        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 3
            assert len(row.result_counts["sample_count"]) == 3
            assert len(row.result_counts["quality_score"]) == 3
            # Check main array: A: 10+5=15, B: 20+0=20, C: 0+15=15
            assert row.result_array[0] == 15
            assert row.result_array[1] == 20
            assert row.result_array[2] == 15
            # Check sample_count: A: 100+50=150, B: 200+0=200, C: 0+150=150
            assert row.result_counts["sample_count"][0] == 150
            assert row.result_counts["sample_count"][1] == 200
            assert row.result_counts["sample_count"][2] == 150
            # Check quality_score: A: 30+25=55, B: 40+0=40, C: 0+35=35
            assert row.result_counts["quality_score"][0] == 55
            assert row.result_counts["quality_score"][1] == 40
            assert row.result_counts["quality_score"][2] == 35

    def test_merge_arrays_with_count_arrays_diff(self, sample_ht):
        """Test merging arrays with count_arrays using diff operation."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
        )

        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="diff",
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        # Expected: only groups from first array (A, B)
        expected_meta = [{"group": "A"}, {"group": "B"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_array) == 2
            assert len(row.result_counts["sample_count"]) == 2
            # Check main array: A: 10-5=5, B: 20-15=5
            assert row.result_array[0] == 5
            assert row.result_array[1] == 5
            # Check count array: A: 100-50=50, B: 200-150=50
            assert row.result_counts["sample_count"][0] == 50
            assert row.result_counts["sample_count"][1] == 50

    def test_merge_arrays_with_count_arrays_validation_errors(self, sample_ht):
        """Test validation errors for count_arrays in merge_array_expressions."""
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
            count1=hl.array([100, 200]),
        )
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Test mismatched count array lengths.
        with pytest.raises(
            ValueError, match="Length of count_array 'test' and meta must be equal!"
        ):
            merge_array_expressions(
                [ht.array1, ht.array2],
                [meta1, meta2],
                count_arrays={"test": [hl.array([1])]},
            )

        # Test count array with different length than main arrays.
        with pytest.raises(
            ValueError,
            match="Length of count_array 'sample_count' and meta must be equal!",
        ):
            merge_array_expressions(
                [ht.array1, ht.array2],
                [meta1, meta2],
                count_arrays={
                    "sample_count": [ht.count1]
                },  # Only one count array instead of two.
            )


class TestMergeFreqArrays:
    """Test the merge_freq_arrays function."""

    @pytest.fixture
    def sample_ht(self):
        """Create a sample Hail Table for testing."""
        return hl.Table.parallelize(
            [
                {"id": 1, "group": "A"},
                {"id": 2, "group": "B"},
                {"id": 3, "group": "C"},
            ],
            hl.tstruct(id=hl.tint32, group=hl.tstr),
        )

    def test_merge_freq_arrays_sum(self, sample_ht):
        """Test merging frequency arrays with sum operation."""
        ht = sample_ht.annotate(
            freq1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            freq2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_freq, result_meta = merge_freq_arrays(
            [ht.freq1, ht.freq2], [meta1, meta2], operation="sum"
        )[:2]

        result = ht.select(result_freq=result_freq).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C with summed values and calculated AF.
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_freq) == 3
            # Check first element (group A): AC=10+5=15, AN=100+50=150, AF=15/150=0.1.
            first_struct = row.result_freq[0]
            assert first_struct.AC == 15
            assert first_struct.AN == 150
            assert first_struct.homozygote_count == 3
            assert first_struct.AF == 0.1

    def test_merge_freq_arrays_diff(self, sample_ht):
        """Test merging frequency arrays with diff operation."""
        ht = sample_ht.annotate(
            freq1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            freq2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_freq, result_meta = merge_freq_arrays(
            [ht.freq1, ht.freq2], [meta1, meta2], operation="diff"
        )[:2]

        result = ht.select(result_freq=result_freq).collect()
        result_meta = hl.eval(result_meta)

        # Expected: only groups from first array (A, B) with subtracted values.
        expected_meta = [{"group": "A"}, {"group": "B"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_freq) == 2
            # Check first element (group A): AC=10-5=5, AN=100-50=50, AF=5/50=0.1.
            first_struct = row.result_freq[0]
            assert first_struct.AC == 5
            assert first_struct.AN == 50
            assert first_struct.homozygote_count == 1
            assert first_struct.AF == 0.1

    def test_merge_freq_arrays_with_count_arrays(self, sample_ht):
        """Test merging frequency arrays with count arrays."""
        ht = sample_ht.annotate(
            freq1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            freq2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        result_freq, result_meta, result_counts = merge_freq_arrays(
            [ht.freq1, ht.freq2],
            [meta1, meta2],
            operation="sum",
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        result = ht.select(
            result_freq=result_freq, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        for row in result:
            assert len(row.result_freq) == 3
            assert len(row.result_counts["sample_count"]) == 3
            # Check count array: A: 100+50=150, B: 200+0=200, C: 0+150=150.
            assert row.result_counts["sample_count"][0] == 150
            assert row.result_counts["sample_count"][1] == 200
            assert row.result_counts["sample_count"][2] == 150

    def test_merge_freq_arrays_with_negatives_set_to_zero(self, sample_ht):
        """Test that negative values in frequency arrays are set to zero."""
        # Create frequency arrays that will result in negative values.
        ht = sample_ht.annotate(
            freq1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            freq2=hl.array(
                [
                    hl.struct(AC=15, AN=50, homozygote_count=3),
                    hl.struct(AC=25, AN=150, homozygote_count=5),
                ]
            ),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        result_freq = merge_freq_arrays(
            [ht.freq1, ht.freq2],
            [meta1, meta2],
            operation="diff",
            set_negatives_to_zero=True,
        )[0]

        result = ht.select(result_freq=result_freq).collect()

        for row in result:
            # A: AC=10-15=-5 -> 0, AN=100-50=50, AF=0/50=0.0.
            first_struct = row.result_freq[0]
            assert first_struct.AC == 0
            assert first_struct.AN == 50
            assert first_struct.AF == 0.0

    def test_merge_freq_arrays_validation_errors(self, sample_ht):
        """Test validation errors in merge_freq_arrays."""
        ht = sample_ht.annotate(
            freq1=hl.array(
                [
                    hl.struct(AC=10, AN=100, homozygote_count=2),
                    hl.struct(AC=20, AN=200, homozygote_count=4),
                ]
            ),
            freq2=hl.array(
                [
                    hl.struct(AC=5, AN=50, homozygote_count=1),
                    hl.struct(AC=15, AN=150, homozygote_count=3),
                ]
            ),
        )

        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        with pytest.raises(
            ValueError, match="Must provide at least two arrays to merge!"
        ):
            merge_freq_arrays([ht.freq1], [meta1])

        with pytest.raises(
            ValueError, match="Length of arrays and meta must be equal!"
        ):
            merge_freq_arrays([ht.freq1, ht.freq2], [meta1])

        with pytest.raises(
            ValueError, match="Operation must be either 'sum' or 'diff'!"
        ):
            merge_freq_arrays([ht.freq1, ht.freq2], [meta1, meta2], operation="invalid")


class TestMergeHistograms:
    """Test the merge_histograms function."""

    @pytest.fixture
    def sample_ht(self):
        """Create a sample Hail Table for testing."""
        return hl.Table.parallelize(
            [
                {"id": 1, "group": "A"},
                {"id": 2, "group": "B"},
                {"id": 3, "group": "C"},
            ],
            hl.tstruct(id=hl.tint32, group=hl.tstr),
        )

    def test_merge_histograms_sum(self, sample_ht):
        """Test merging histograms with sum operation."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([3, 7, 12]),
                n_smaller=1,
                n_larger=5,
            ),
        )

        result_hist = merge_histograms([ht.hist1, ht.hist2], operation="sum")
        result = ht.select(result_hist=result_hist).collect()

        for row in result:
            assert row.result_hist.bin_edges == [0, 10, 20, 30]
            # Check bin_freq (summed): [5+3, 10+7, 15+12] = [8, 17, 27].
            assert row.result_hist.bin_freq == [8, 17, 27]
            assert row.result_hist.n_smaller == 3
            assert row.result_hist.n_larger == 13

    def test_merge_histograms_diff(self, sample_ht):
        """Test merging histograms with diff operation."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([10, 20, 30]),
                n_smaller=5,
                n_larger=15,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([3, 7, 12]),
                n_smaller=2,
                n_larger=8,
            ),
        )

        result_hist = merge_histograms([ht.hist1, ht.hist2], operation="diff")
        result = ht.select(result_hist=result_hist).collect()

        for row in result:
            assert row.result_hist.bin_edges == [0, 10, 20, 30]
            # Check bin_freq (subtracted): [10-3, 20-7, 30-12] = [7, 13, 18].
            assert row.result_hist.bin_freq == [7, 13, 18]
            assert row.result_hist.n_smaller == 3
            assert row.result_hist.n_larger == 7

    def test_merge_histograms_with_negatives_error(self, sample_ht):
        """Test that negative values raise error when set_negatives_to_zero=False."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([8, 15, 20]),
                n_smaller=5,
                n_larger=12,
            ),
        )

        # Should raise an error due to negative values in diff operation.
        result_hist = merge_histograms([ht.hist1, ht.hist2], operation="diff")
        with pytest.raises(Exception):
            ht.select(result_hist=result_hist).collect()

    def test_merge_histograms_with_negatives_set_to_zero(self, sample_ht):
        """Test that negative values are set to zero when set_negatives_to_zero=True."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([8, 15, 20]),
                n_smaller=5,
                n_larger=12,
            ),
        )

        result_hist = merge_histograms(
            [ht.hist1, ht.hist2], operation="diff", set_negatives_to_zero=True
        )
        result = ht.select(result_hist=result_hist).collect()

        for row in result:
            # bin_freq: [5-8, 10-15, 15-20] = [-3, -5, -5] -> [0, 0, 0].
            assert row.result_hist.bin_freq == [0, 0, 0]
            # n_smaller: 2 - 5 = -3 -> 0.
            assert row.result_hist.n_smaller == 0
            # n_larger: 8 - 12 = -4 -> 0.
            assert row.result_hist.n_larger == 0

    def test_merge_histograms_validation_errors(self, sample_ht):
        """Test that merge_histograms raises appropriate ValueError for invalid inputs."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([3, 7, 12]),
                n_smaller=1,
                n_larger=5,
            ),
        )

        with pytest.raises(
            ValueError, match="Must provide at least two histograms to merge!"
        ):
            merge_histograms([ht.hist1])

        with pytest.raises(
            ValueError, match="Operation must be either 'sum' or 'diff'!"
        ):
            merge_histograms([ht.hist1, ht.hist2], operation="invalid")

    def test_merge_histograms_multiple_histograms(self, sample_ht):
        """Test merging multiple histograms."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([3, 7, 12]),
                n_smaller=1,
                n_larger=5,
            ),
            hist3=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([2, 4, 6]),
                n_smaller=0,
                n_larger=3,
            ),
        )

        result_hist = merge_histograms([ht.hist1, ht.hist2, ht.hist3], operation="sum")
        result = ht.select(result_hist=result_hist).collect()

        for row in result:
            assert row.result_hist.bin_edges == [0, 10, 20, 30]
            # Check bin_freq (summed): [5+3+2, 10+7+4, 15+12+6] = [10, 21, 33].
            assert row.result_hist.bin_freq == [10, 21, 33]
            assert row.result_hist.n_smaller == 3
            assert row.result_hist.n_larger == 16

    def test_merge_histograms_with_missing_values(self, sample_ht):
        """Test merging histograms with missing values."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([3, None, 12]),
                n_smaller=1,
                n_larger=5,
            ),
        )

        result_hist = merge_histograms([ht.hist1, ht.hist2], operation="sum")
        result = ht.select(result_hist=result_hist).collect()

        for row in result:
            # Check bin_freq: [5+3, 10+0, 15+12] = [8, 10, 27] (missing treated as 0).
            assert row.result_hist.bin_freq == [8, 10, 27]
            assert row.result_hist.n_smaller == 3
            assert row.result_hist.n_larger == 13

    def test_merge_histograms_sum_with_negatives_set_to_zero(self, sample_ht):
        """Test that negative values are handled correctly in sum operation when set_negatives_to_zero=True."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([-3, -5, -7]),
                n_smaller=-1,
                n_larger=-3,
            ),
        )

        # Test that negative values are set to zero when set_negatives_to_zero=True
        result_hist = merge_histograms(
            [ht.hist1, ht.hist2], operation="sum", set_negatives_to_zero=True
        )
        result = ht.select(result_hist=result_hist).collect()

        for row in result:
            # bin_freq: [5+(-3), 10+(-5), 15+(-7)] = [2, 5, 8] (no negatives in sum)
            assert row.result_hist.bin_freq == [2, 5, 8]
            # n_smaller: 2 + (-1) = 1
            assert row.result_hist.n_smaller == 1
            # n_larger: 8 + (-3) = 5
            assert row.result_hist.n_larger == 5

    def test_merge_histograms_sum_with_negatives_error(self, sample_ht):
        """Test that negative values in sum operation raise error when set_negatives_to_zero=False."""
        ht = sample_ht.annotate(
            hist1=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([5, 10, 15]),
                n_smaller=2,
                n_larger=8,
            ),
            hist2=hl.struct(
                bin_edges=hl.array([0, 10, 20, 30]),
                bin_freq=hl.array([-8, -15, -20]),
                n_smaller=-5,
                n_larger=-12,
            ),
        )

        # Test that negative values in sum operation raise error when
        # set_negatives_to_zero=False
        result_hist = merge_histograms([ht.hist1, ht.hist2], operation="sum")
        with pytest.raises(Exception):
            ht.select(result_hist=result_hist).collect()


class TestVRSFunctions:
    """Test the VRS-related functions."""

    # Input values are gnomAD exomes v4.1 GRCh38 chr16 VRS-annotated variants.
    # Generated with vrs-python (ga4gh.vrs==2.2.0)
    CHR16_SNP = {
        "contig": "chr16",
        "pos": 30223487,
        "ref": "G",
        "alt": "C",
        # VCF INFO: VRS_Allele_IDs=ga4gh:VA.8pQgh_...,ga4gh:VA.Ljkc...
        "vrs_ref_allele_id": "ga4gh:VA.8pQgh_5I1jlQdH5SVfXcZnQv-QT3wq2m",
        "vrs_alt_allele_id": "ga4gh:VA.Ljkcleg2iX7DazIE3lDKIx2tAuuIJilw",
        # VCF INFO: VRS_Starts=30223486,30223486; VRS_Ends=30223487,30223487
        "vrs_ref_start": 30223486,
        "vrs_alt_start": 30223486,
        "vrs_ref_end": 30223487,
        "vrs_alt_end": 30223487,
    }
    CHR16_DELETION_RLE = {
        "contig": "chr16",
        "pos": 30223488,
        "ref": "AGCTG",
        "alt": "A",
        # VCF INFO: VRS_Allele_IDs=ga4gh:VA.H-Kn...,ga4gh:VA.o0wp...
        "vrs_ref_allele_id": "ga4gh:VA.H-KnDD4KLqQ8hAC1lHo2jur_NHFKJV6I",
        "vrs_alt_allele_id": "ga4gh:VA.o0wpSwOwTnYfNz3Yz1E4AzZGQYyBg6e3",
        # VCF INFO: VRS_Starts=30223487,30223488; VRS_Ends=30223492,30223492
        "vrs_ref_start": 30223487,
        "vrs_alt_start": 30223488,
        "vrs_ref_end": 30223492,
        "vrs_alt_end": 30223492,
        "rle_length": 0,
        "rle_repeat_subunit_length": 4,
    }
    CHR16_RLE_NOT_DELETION = {
        "contig": "chr16",
        "pos": 13828,
        "ref": "T",
        "alt": "TA",
        # VCF INFO: VRS_Allele_IDs=ga4gh:VA.GmRX-...,ga4gh:VA._k6Xx...
        "vrs_ref_allele_id": "ga4gh:VA.GmRX-v5NmDfEaY24aOPnhDUQoZwJszKZ",
        "vrs_alt_allele_id": "ga4gh:VA._k6XxTZDGXxhe0I-OVk6Tq8AOdVa-QYN",
        # VCF INFO: VRS_Starts=13827,13828; VRS_Ends=13828,13829
        "vrs_ref_start": 13827,
        "vrs_alt_start": 13828,
        "vrs_ref_end": 13828,
        "vrs_alt_end": 13829,
        "rle_length": 2,
        "rle_repeat_subunit_length": 1,
        "rle_sequence": "AA",
    }
    CHR16_RLE_REPEAT_EXPANSION = {
        "contig": "chr16",
        "pos": 12247,
        "ref": "C",
        "alt": "CTGG",
        # VCF INFO: VRS_Allele_IDs=ga4gh:VA.dFms...,ga4gh:VA.yQ4C...
        "vrs_ref_allele_id": "ga4gh:VA.dFmsD4FFHwdOJbPV_6armbHSKs9FWQec",
        "vrs_alt_allele_id": "ga4gh:VA.yQ4CuV2WkXfUed8UU8dwzf6POVL_VXPg",
        # VCF INFO: VRS_Starts=12246,12247; VRS_Ends=12247,12252
        "vrs_ref_start": 12246,
        "vrs_alt_start": 12247,
        "vrs_ref_end": 12247,
        "vrs_alt_end": 12252,
        # VCF INFO: VRS_States=C,TGGTGGTG; VRS_Lengths=1,8; VRS_RepeatSubunitLengths=1,3
        "rle_length": 8,
        "rle_repeat_subunit_length": 3,
        "rle_sequence": "TGGTGGTG",
    }

    def test_vrs_identifier_generation(self):
        """Test that GA4GH identifiers are generated correctly."""
        # Test that we can import and access VRS modules
        assert hasattr(ga4gh_core, "__version__")
        assert hasattr(ga4gh_vrs, "models")
        assert hasattr(ga4gh_vrs.models, "SequenceLocation")

        # Test that VRS 2.0.1+ API is available
        assert hasattr(
            ga4gh_core, "ga4gh_identify"
        ), "VRS 2.0.1+ ga4gh_identify function not found"

        # Test that we can create a VRS object using the 2.0.1+ API
        seq_loc = ga4gh_vrs.models.SequenceLocation(
            sequenceReference="ga4gh:SQ.test",
            start=1,
            end=2,
        )
        assert seq_loc is not None

    def test_vrs_error_handling(self):
        """Test that VRS functions handle errors gracefully."""
        # Test with invalid locus (should raise appropriate error)
        with pytest.raises((AttributeError, TypeError)):
            # This should fail because we're not providing proper Hail locus objects
            add_gks_vrs("invalid_locus", "invalid_vrs")

    def test_add_gks_vrs_actual_api_call(self):
        """Test that add_gks_vrs calls the VRS API with real SNP data."""
        locus = hl.locus(
            self.CHR16_SNP["contig"],
            self.CHR16_SNP["pos"],
            reference_genome="GRCh38",
        )
        vrs_struct = hl.struct(
            VRS_Allele_IDs=[
                self.CHR16_SNP["vrs_ref_allele_id"],
                self.CHR16_SNP["vrs_alt_allele_id"],
            ],
            VRS_Starts=[
                self.CHR16_SNP["vrs_ref_start"],
                self.CHR16_SNP["vrs_alt_start"],
            ],
            VRS_Ends=[self.CHR16_SNP["vrs_ref_end"], self.CHR16_SNP["vrs_alt_end"]],
            VRS_States=[self.CHR16_SNP["ref"], self.CHR16_SNP["alt"]],
            VRS_Lengths=hl.literal([1, None], hl.tarray(hl.tint32)),
            VRS_RepeatSubunitLengths=hl.literal([1, None], hl.tarray(hl.tint32)),
        )

        # Evaluate to get actual Python objects
        locus_py = hl.eval(locus)
        vrs_py = hl.eval(vrs_struct)

        result = add_gks_vrs(locus_py, vrs_py)
        # Validate with Pydantic, but do not renormalize.
        assert isinstance(result, list)
        assert len(result) == 2
        _ = ga4gh_vrs.models.Allele(**result[0])
        _ = ga4gh_vrs.models.Allele(**result[1])
        alt_result = result[1]

        # Verify the result has the expected VRS structure
        assert alt_result["type"] == "Allele"
        assert "location" in alt_result
        assert alt_result["location"]["type"] == "SequenceLocation"
        assert "id" in alt_result["location"]  # This comes from ga4gh_identify() call
        assert "state" in alt_result
        assert alt_result["state"]["type"] == "LiteralSequenceExpression"
        assert alt_result["state"]["sequence"] == self.CHR16_SNP["alt"]

        # Verify the chromosome ID matches our VRS_CHROM_IDS mapping
        expected_chr16_id = VRS_CHROM_IDS["GRCh38"]["chr16"]
        assert expected_chr16_id == "SQ.yC_0RBj3fgBlvgyAuycbzdubtLxq-rE0"
        assert (
            alt_result["location"]["sequenceReference"]["refgetAccession"]
            == expected_chr16_id
        )

        # Verify the location ID was generated by the VRS API
        location_id = alt_result["location"]["id"]
        assert isinstance(location_id, str)
        assert location_id.startswith(
            "ga4gh:SL"
        )  # VRS SequenceLocation identifiers start with ga4gh:SL

    def test_add_gks_vrs_reference_length_expression_deletion(self):
        """Test that add_gks_vrs returns a valid ReferenceLengthExpression for a deletion."""
        locus = hl.locus(
            self.CHR16_DELETION_RLE["contig"],
            self.CHR16_DELETION_RLE["pos"],
            reference_genome="GRCh38",
        )

        vrs_struct = hl.struct(
            VRS_Allele_IDs=[
                self.CHR16_DELETION_RLE["vrs_ref_allele_id"],
                self.CHR16_DELETION_RLE["vrs_alt_allele_id"],
            ],
            VRS_Starts=[
                self.CHR16_DELETION_RLE["vrs_ref_start"],
                self.CHR16_DELETION_RLE["vrs_alt_start"],
            ],
            VRS_Ends=[
                self.CHR16_DELETION_RLE["vrs_ref_end"],
                self.CHR16_DELETION_RLE["vrs_alt_end"],
            ],
            VRS_States=hl.literal(
                [self.CHR16_DELETION_RLE["ref"], None], hl.tarray(hl.tstr)
            ),
            VRS_Lengths=hl.literal([5, 0], hl.tarray(hl.tint32)),
            VRS_RepeatSubunitLengths=hl.literal([5, 4], hl.tarray(hl.tint32)),
        )

        locus_py = hl.eval(locus)
        vrs_py = hl.eval(vrs_struct)
        result = add_gks_vrs(locus_py, vrs_py)

        assert isinstance(result, list)
        assert len(result) == 2
        alt_result = result[1]
        _ = ga4gh_vrs.models.Allele(**alt_result)
        assert alt_result["state"]["type"] == "ReferenceLengthExpression"
        assert alt_result["state"]["length"] == self.CHR16_DELETION_RLE["rle_length"]
        assert (
            alt_result["state"]["repeatSubunitLength"]
            == self.CHR16_DELETION_RLE["rle_repeat_subunit_length"]
        )

    def test_add_gks_vrs_reference_length_expression_not_deletion(self):
        """Test a ReferenceLengthExpression where the alternate length is non-zero."""
        locus = hl.locus(
            self.CHR16_RLE_NOT_DELETION["contig"],
            self.CHR16_RLE_NOT_DELETION["pos"],
            reference_genome="GRCh38",
        )

        vrs_struct = hl.struct(
            VRS_Allele_IDs=[
                self.CHR16_RLE_NOT_DELETION["vrs_ref_allele_id"],
                self.CHR16_RLE_NOT_DELETION["vrs_alt_allele_id"],
            ],
            VRS_Starts=[
                self.CHR16_RLE_NOT_DELETION["vrs_ref_start"],
                self.CHR16_RLE_NOT_DELETION["vrs_alt_start"],
            ],
            VRS_Ends=[
                self.CHR16_RLE_NOT_DELETION["vrs_ref_end"],
                self.CHR16_RLE_NOT_DELETION["vrs_alt_end"],
            ],
            VRS_States=[
                self.CHR16_RLE_NOT_DELETION["ref"],
                self.CHR16_RLE_NOT_DELETION["rle_sequence"],
            ],
            VRS_Lengths=hl.literal([1, 2], hl.tarray(hl.tint32)),
            VRS_RepeatSubunitLengths=hl.literal([1, 1], hl.tarray(hl.tint32)),
        )

        locus_py = hl.eval(locus)
        vrs_py = hl.eval(vrs_struct)
        result = add_gks_vrs(locus_py, vrs_py)

        assert isinstance(result, list)
        assert len(result) == 2
        alt_result = result[1]
        _ = ga4gh_vrs.models.Allele(**alt_result)
        assert alt_result["state"]["type"] == "ReferenceLengthExpression"
        assert (
            alt_result["state"]["length"] == self.CHR16_RLE_NOT_DELETION["rle_length"]
        )
        assert (
            alt_result["state"]["repeatSubunitLength"]
            == self.CHR16_RLE_NOT_DELETION["rle_repeat_subunit_length"]
        )
        assert (
            alt_result["state"]["sequence"]
            == self.CHR16_RLE_NOT_DELETION["rle_sequence"]
        )

    def test_add_gks_vrs_reference_length_expression_repeat_expansion(self):
        """Test a ReferenceLengthExpression where the ALT is longer than the REF."""
        assert len(self.CHR16_RLE_REPEAT_EXPANSION["alt"]) > len(
            self.CHR16_RLE_REPEAT_EXPANSION["ref"]
        )

        locus = hl.locus(
            self.CHR16_RLE_REPEAT_EXPANSION["contig"],
            self.CHR16_RLE_REPEAT_EXPANSION["pos"],
            reference_genome="GRCh38",
        )

        vrs_struct = hl.struct(
            VRS_Allele_IDs=[
                self.CHR16_RLE_REPEAT_EXPANSION["vrs_ref_allele_id"],
                self.CHR16_RLE_REPEAT_EXPANSION["vrs_alt_allele_id"],
            ],
            VRS_Starts=[
                self.CHR16_RLE_REPEAT_EXPANSION["vrs_ref_start"],
                self.CHR16_RLE_REPEAT_EXPANSION["vrs_alt_start"],
            ],
            VRS_Ends=[
                self.CHR16_RLE_REPEAT_EXPANSION["vrs_ref_end"],
                self.CHR16_RLE_REPEAT_EXPANSION["vrs_alt_end"],
            ],
            VRS_States=[
                self.CHR16_RLE_REPEAT_EXPANSION["ref"],
                self.CHR16_RLE_REPEAT_EXPANSION["rle_sequence"],
            ],
            VRS_Lengths=hl.literal(
                [1, self.CHR16_RLE_REPEAT_EXPANSION["rle_length"]], hl.tarray(hl.tint32)
            ),
            VRS_RepeatSubunitLengths=hl.literal(
                [1, self.CHR16_RLE_REPEAT_EXPANSION["rle_repeat_subunit_length"]],
                hl.tarray(hl.tint32),
            ),
        )

        locus_py = hl.eval(locus)
        vrs_py = hl.eval(vrs_struct)
        result = add_gks_vrs(locus_py, vrs_py)

        assert isinstance(result, list)
        assert len(result) == 2
        alt_result = result[1]
        _ = ga4gh_vrs.models.Allele(**alt_result)
        assert alt_result["state"]["type"] == "ReferenceLengthExpression"
        assert (
            alt_result["state"]["length"]
            == self.CHR16_RLE_REPEAT_EXPANSION["rle_length"]
        )
        assert (
            alt_result["state"]["repeatSubunitLength"]
            == self.CHR16_RLE_REPEAT_EXPANSION["rle_repeat_subunit_length"]
        )
        assert (
            alt_result["state"]["sequence"]
            == self.CHR16_RLE_REPEAT_EXPANSION["rle_sequence"]
        )

    def test_add_gks_vrs_grch37_chromosome_mapping(self):
        """Test that VRS chromosome mapping works correctly for GRCh37."""
        # Test GRCh37 (note different chromosome naming: "1" vs "chr1")
        locus_grch37 = hl.locus("1", 100, reference_genome="GRCh37")
        vrs_struct = hl.struct(
            VRS_Allele_IDs=["ga4gh:VA." + "A" * 32, "ga4gh:VA." + "B" * 32],
            VRS_Starts=[99, 99],
            VRS_Ends=[100, 100],
            VRS_States=["A", "T"],
            VRS_Lengths=[1, 1],
            VRS_RepeatSubunitLengths=hl.literal([None, None], hl.tarray(hl.tint32)),
        )

        locus_py = hl.eval(locus_grch37)
        vrs_py = hl.eval(vrs_struct)

        result = add_gks_vrs(locus_py, vrs_py)
        # Validate with Pydantic, but do not renormalize.
        assert isinstance(result, list)
        assert len(result) == 2
        _ = ga4gh_vrs.models.Allele(**result[0])
        _ = ga4gh_vrs.models.Allele(**result[1])

        # Verify GRCh37 chromosome 1 mapping
        expected_chr1_grch37_id = VRS_CHROM_IDS["GRCh37"]["1"]
        assert expected_chr1_grch37_id == "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU"
        assert (
            result[1]["location"]["sequenceReference"]["refgetAccession"]
            == expected_chr1_grch37_id
        )

        # Verify the API was called and generated a proper identifier
        location_id = result[1]["location"]["id"]
        assert isinstance(location_id, str)
        assert location_id.startswith(
            "ga4gh:SL"
        )  # VRS SequenceLocation identifiers start with ga4gh:SL


class TestAnnotateDownsamplings:
    """Test the annotate_downsamplings function."""

    @pytest.fixture
    def sample_matrix_table(self):
        """Create a sample MatrixTable for testing."""
        samples = [
            {"s": "sample1", "gen_anc": "AFR"},
            {"s": "sample2", "gen_anc": "EUR"},
            {"s": "sample3", "gen_anc": "AFR"},
            {"s": "sample4", "gen_anc": "EUR"},
            {"s": "sample5", "gen_anc": "SAS"},
        ]

        variants = [
            {
                "locus": hl.locus("chr1", 1000, reference_genome="GRCh38"),
                "alleles": ["A", "T"],
            },
            {
                "locus": hl.locus("chr1", 2000, reference_genome="GRCh38"),
                "alleles": ["C", "G"],
            },
            {
                "locus": hl.locus("chr1", 3000, reference_genome="GRCh38"),
                "alleles": ["T", "A"],
            },
        ]

        sample_table = hl.Table.parallelize(
            samples,
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr),
        ).key_by("s")

        entries = []
        for variant in variants:
            for sample in samples:
                entries.append(
                    {
                        "locus": variant["locus"],
                        "alleles": variant["alleles"],
                        "s": sample["s"],
                        "GT": hl.call(0, 1),
                    }
                )

        mt = hl.Table.parallelize(
            entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])

        mt = mt.annotate_cols(gen_anc=sample_table[mt.s].gen_anc)

        return mt

    @pytest.fixture
    def sample_table(self):
        """Create a sample Table for testing."""
        return hl.Table.parallelize(
            [
                {"s": "sample1", "sex": "XX", "gen_anc": "AFR"},
                {"s": "sample2", "sex": "XY", "gen_anc": "EUR"},
                {"s": "sample3", "sex": "XX", "gen_anc": "AFR"},
                {"s": "sample4", "sex": "XY", "gen_anc": "EUR"},
                {"s": "sample5", "sex": "XX", "gen_anc": "SAS"},
            ],
            hl.tstruct(s=hl.tstr, sex=hl.tstr, gen_anc=hl.tstr),
        ).key_by("s")

    def test_annotate_downsamplings_matrix_table_no_gen_anc(self, sample_matrix_table):
        """Test annotate_downsamplings with MatrixTable input without genetic ancestry."""
        downsamplings = [2, 3, 4]

        result = annotate_downsamplings(sample_matrix_table, downsamplings)

        assert isinstance(result, hl.MatrixTable)
        assert "downsampling" in result.col.dtype
        assert "downsamplings" in result.globals.dtype

        result_downsamplings = hl.eval(result.downsamplings)
        assert result_downsamplings == [2, 3, 4]

        sample_cols = result.cols().collect()
        for col in sample_cols:
            assert "global_idx" in col.downsampling

    def test_annotate_downsamplings_matrix_table_with_gen_anc(
        self, sample_matrix_table
    ):
        """Test annotate_downsamplings with MatrixTable input with genetic ancestry."""
        downsamplings = [2, 3, 4]
        gen_anc_expr = sample_matrix_table.gen_anc

        result = annotate_downsamplings(
            sample_matrix_table, downsamplings, gen_anc_expr
        )

        assert isinstance(result, hl.MatrixTable)
        assert "downsampling" in result.col.dtype
        assert "downsamplings" in result.globals.dtype
        assert "ds_gen_anc_counts" in result.globals.dtype

        result_downsamplings = hl.eval(result.downsamplings)
        assert result_downsamplings == [1, 2, 3, 4]

        gen_anc_counts = hl.eval(result.ds_gen_anc_counts)
        assert gen_anc_counts == {"AFR": 2, "EUR": 2, "SAS": 1}

        sample_cols = result.cols().collect()
        for col in sample_cols:
            assert "global_idx" in col.downsampling
            assert "gen_anc_idx" in col.downsampling

    def test_annotate_downsamplings_table_no_gen_anc(self, sample_table):
        """Test annotate_downsamplings with Table input without genetic ancestry."""
        downsamplings = [2, 3, 4]

        result = annotate_downsamplings(sample_table, downsamplings)

        assert isinstance(result, hl.Table)
        assert "downsampling" in result.row.dtype
        assert "downsamplings" in result.globals.dtype

        result_downsamplings = hl.eval(result.downsamplings)
        assert result_downsamplings == [2, 3, 4]

        rows = result.collect()
        for row in rows:
            assert "global_idx" in row.downsampling

    def test_annotate_downsamplings_table_with_gen_anc(self, sample_table):
        """Test annotate_downsamplings with Table input with genetic ancestry."""
        downsamplings = [2, 3, 4]
        gen_anc_expr = sample_table.gen_anc

        result = annotate_downsamplings(sample_table, downsamplings, gen_anc_expr)

        assert isinstance(result, hl.Table)
        assert "downsampling" in result.row.dtype
        assert "downsamplings" in result.globals.dtype
        assert "ds_gen_anc_counts" in result.globals.dtype

        result_downsamplings = hl.eval(result.downsamplings)
        assert result_downsamplings == [1, 2, 3, 4]

        gen_anc_counts = hl.eval(result.ds_gen_anc_counts)
        assert gen_anc_counts == {"AFR": 2, "EUR": 2, "SAS": 1}

        rows = result.collect()
        for row in rows:
            assert "global_idx" in row.downsampling
            assert "gen_anc_idx" in row.downsampling

    @pytest.fixture
    def gen_anc_sized_table(self):
        """Create a Table with distinct gen-anc group sizes (AFR=5, EUR=4, SAS=3)."""
        rows = (
            [{"s": f"a{i}", "gen_anc": "AFR"} for i in range(5)]
            + [{"s": f"e{i}", "gen_anc": "EUR"} for i in range(4)]
            + [{"s": f"x{i}", "gen_anc": "SAS"} for i in range(3)]
        )
        return hl.Table.parallelize(
            rows, hl.tstruct(s=hl.tstr, gen_anc=hl.tstr)
        ).key_by("s")

    def test_gen_anc_group_sizes_added_dynamically(self, gen_anc_sized_table):
        """Test that per-group sizes are appended to the downsamplings list.

        The input list holds only a "global" size (6); each genetic ancestry
        group's sample count (AFR=5, EUR=4, SAS=3) is added dynamically by passing
        `gen_anc_expr` to the function so per-group sizes do not need to be included
        in the DOWNSAMPLINGS constant.
        """
        result = annotate_downsamplings(
            gen_anc_sized_table, [6], gen_anc_expr=gen_anc_sized_table.gen_anc
        )

        # None of the group sizes (3, 4, 5) were in the input list yet all appear.
        assert hl.eval(result.downsamplings) == [3, 4, 5, 6]
        assert hl.eval(result.ds_gen_anc_counts) == {"AFR": 5, "EUR": 4, "SAS": 3}

    def test_gen_ancs_to_downsample_restricts_added_groups(self, gen_anc_sized_table):
        """Test that `gen_ancs_to_downsample` limits which group sizes are added."""
        result = annotate_downsamplings(
            gen_anc_sized_table,
            [6],
            gen_anc_expr=gen_anc_sized_table.gen_anc,
            gen_ancs_to_downsample=["AFR", "SAS"],
        )

        # EUR is excluded, so its size (4) is neither added to the list nor
        # counted; AFR (5) and SAS (3) are.
        assert hl.eval(result.downsamplings) == [3, 5, 6]
        assert hl.eval(result.ds_gen_anc_counts) == {"AFR": 5, "SAS": 3}

    def test_downsamplings_exceeding_total_are_dropped(self, gen_anc_sized_table):
        """Test that requested sizes larger than the sample total are dropped."""
        result = annotate_downsamplings(
            gen_anc_sized_table, [6, 100], gen_anc_expr=gen_anc_sized_table.gen_anc
        )

        # 100 > 12 total samples, so it is dropped; the dynamic group sizes remain.
        assert hl.eval(result.downsamplings) == [3, 4, 5, 6]

    def test_global_downsampling_kept_when_no_group_reaches_it(
        self, gen_anc_sized_table
    ):
        """Test that a global downsampling is filtered against the full dataset.

        A global downsampling larger than every individual genetic ancestry group
        (and larger than the restricted `gen_ancs_to_downsample` subset total) is
        retained as long as it does not exceed the full dataset. Guards the
        `total_gen_anc_counts` filtering fix.
        """
        # Groups AFR=5, EUR=4, SAS=3 (full total 12). 9 exceeds every group and
        # the AFR+SAS subset total (8), but is <= 12, so it must be kept.
        result = annotate_downsamplings(
            gen_anc_sized_table,
            [9],
            gen_anc_expr=gen_anc_sized_table.gen_anc,
            gen_ancs_to_downsample=["AFR", "SAS"],
        )

        assert hl.eval(result.downsamplings) == [3, 5, 9]

    def test_gen_ancs_to_downsample_unknown_group_raises(self, gen_anc_sized_table):
        """Test that requesting a group absent from the data raises a ValueError."""
        with pytest.raises(ValueError, match="not present in the data"):
            annotate_downsamplings(
                gen_anc_sized_table,
                [6],
                gen_anc_expr=gen_anc_sized_table.gen_anc,
                gen_ancs_to_downsample=["AFR", "TYPO"],
            )

    def test_gen_ancs_to_downsample_requires_gen_anc_expr(self, gen_anc_sized_table):
        """Test that `gen_ancs_to_downsample` without `gen_anc_expr` raises."""
        with pytest.raises(ValueError, match="requires `gen_anc_expr`"):
            annotate_downsamplings(
                gen_anc_sized_table, [6], gen_ancs_to_downsample=["AFR"]
            )

    def test_global_idx_is_unique_ranking(self, gen_anc_sized_table):
        """Test that `global_idx` assigns each sample a unique rank 0..n-1."""
        result = annotate_downsamplings(gen_anc_sized_table, [2, 3])

        global_idxs = sorted(r.downsampling.global_idx for r in result.collect())
        # 12 samples get the ranks 0..11 exactly once each.
        assert global_idxs == list(range(gen_anc_sized_table.count()))

    def test_gen_anc_idx_is_per_group_ranking(self, gen_anc_sized_table):
        """Test that `gen_anc_idx` ranks samples 0..(group_size-1) within a group."""
        result = annotate_downsamplings(
            gen_anc_sized_table, [2], gen_anc_expr=gen_anc_sized_table.gen_anc
        )

        per_group = {}
        for r in result.collect():
            per_group.setdefault(r.gen_anc, []).append(r.downsampling.gen_anc_idx)

        # Each group's indices are a contiguous 0..(size-1) ranking.
        assert {g: sorted(v) for g, v in per_group.items()} == {
            "AFR": [0, 1, 2, 3, 4],
            "EUR": [0, 1, 2, 3],
            "SAS": [0, 1, 2],
        }

    def test_downsamplings_deduplicated_and_sorted_with_gen_anc(
        self, gen_anc_sized_table
    ):
        """Test that duplicate/unsorted downsamplings are normalized with gen_anc."""
        result = annotate_downsamplings(
            gen_anc_sized_table,
            [3, 2, 2, 3],
            gen_anc_expr=gen_anc_sized_table.gen_anc,
        )

        # Input deduped to {2, 3} and merged with the sorted group sizes {3, 4, 5}.
        assert hl.eval(result.downsamplings) == [2, 3, 4, 5]

    def test_downsamplings_not_normalized_without_gen_anc(self, gen_anc_sized_table):
        """Document that without `gen_anc_expr` the downsamplings list is used as-is.

        Unlike the `gen_anc_expr` path (which dedupes, sorts, and drops sizes
        larger than the dataset), the no-gen_anc path returns `downsamplings`
        verbatim. This is a characterization test of current behavior.
        """
        # Duplicates and unsorted order are preserved, and 1000 (> 12 samples) is
        # not dropped.
        result = annotate_downsamplings(gen_anc_sized_table, [3, 2, 2, 1000])
        assert hl.eval(result.downsamplings) == [3, 2, 2, 1000]

    def test_downsampling_struct_fields(self, gen_anc_sized_table):
        """Test the per-sample downsampling struct fields and completeness."""
        no_gen_anc = annotate_downsamplings(gen_anc_sized_table, [2, 3])
        assert set(no_gen_anc.downsampling.dtype) == {"global_idx"}

        with_gen_anc = annotate_downsamplings(
            gen_anc_sized_table, [2, 3], gen_anc_expr=gen_anc_sized_table.gen_anc
        )
        assert set(with_gen_anc.downsampling.dtype) == {"global_idx", "gen_anc_idx"}
        # Every sample is annotated with both indices (none missing).
        assert all(
            r.downsampling.global_idx is not None
            and r.downsampling.gen_anc_idx is not None
            for r in with_gen_anc.collect()
        )


class TestGenerateFreqGroupMembershipArray:
    """Test the generate_freq_group_membership_array function."""

    def test_basic_strata_group_membership(self):
        """Test freq_meta and group_membership shape for simple gen-anc strata."""
        rows = [{"s": f"a{i}", "gen_anc": "AFR"} for i in range(3)] + [
            {"s": f"e{i}", "gen_anc": "EUR"} for i in range(2)
        ]
        ht = hl.Table.parallelize(rows, hl.tstruct(s=hl.tstr, gen_anc=hl.tstr)).key_by(
            "s"
        )

        result = generate_freq_group_membership_array(ht, [{"gen_anc": ht.gen_anc}])
        freq_meta = hl.eval(result.freq_meta)

        # Both genetic ancestry groups get a stratum.
        assert {m["gen_anc"] for m in freq_meta if "gen_anc" in m} == {"AFR", "EUR"}
        # group_membership has one boolean per freq_meta entry.
        assert all(len(r.group_membership) == len(freq_meta) for r in result.collect())

    @pytest.fixture
    def downsampling_ht(self):
        """Return a Table with per-group downsamplings (AFR=5, EUR=4, SAS=3).

        Built via `annotate_downsamplings` restricted to AFR and SAS so that EUR
        is present in the data but absent from `ds_gen_anc_counts`.
        """
        rows = (
            [{"s": f"a{i}", "gen_anc": "AFR"} for i in range(5)]
            + [{"s": f"e{i}", "gen_anc": "EUR"} for i in range(4)]
            + [{"s": f"x{i}", "gen_anc": "SAS"} for i in range(3)]
        )
        ht = hl.Table.parallelize(rows, hl.tstruct(s=hl.tstr, gen_anc=hl.tstr)).key_by(
            "s"
        )
        return annotate_downsamplings(
            ht, [6], gen_anc_expr=ht.gen_anc, gen_ancs_to_downsample=["AFR", "SAS"]
        )

    def test_per_group_downsampling_skipped_when_excluded_or_oversized(
        self, downsampling_ht
    ):
        """Test that per-group downsampling strata are skipped appropriately.

        A per-group downsampling stratum is created only when the group is in
        `ds_gen_anc_counts` and the size does not exceed the group's sample count.
        With downsamplings [3, 5, 6] and counts {AFR: 5, SAS: 3}:
          - EUR is excluded entirely (absent from `ds_gen_anc_counts`),
          - AFR keeps 3 and 5 but skips 6 (> 5),
          - SAS keeps only 3 (5 and 6 > 3).
        """
        downsamplings = hl.eval(downsampling_ht.downsamplings)
        ds_gen_anc_counts = hl.eval(downsampling_ht.ds_gen_anc_counts)
        strata_expr = [
            {
                "gen_anc": downsampling_ht.gen_anc,
                "downsampling": downsampling_ht.downsampling,
            }
        ]

        result = generate_freq_group_membership_array(
            downsampling_ht,
            strata_expr,
            downsamplings=downsamplings,
            ds_gen_anc_counts=ds_gen_anc_counts,
        )
        freq_meta = hl.eval(result.freq_meta)

        # Per-group (non-"global") downsampling strata that were created.
        per_group = sorted(
            (m["gen_anc"], int(m["downsampling"]))
            for m in freq_meta
            if "downsampling" in m and m.get("gen_anc") not in (None, "global")
        )

        assert per_group == [("AFR", 3), ("AFR", 5), ("SAS", 3)]
        # EUR is absent from ds_gen_anc_counts, so it gets no per-group stratum.
        assert not any(m.get("gen_anc") == "EUR" for m in freq_meta)

    @pytest.fixture
    def strata_table(self):
        """Return a 3-sample Table (AFR/XX, AFR/XY, EUR/XX) with gen_anc and sex.

        EUR/XY is intentionally absent so a gen_anc x sex cross-product yields a
        zero-sample group.
        """
        return hl.Table.parallelize(
            [
                {"s": "s1", "gen_anc": "AFR", "sex": "XX"},
                {"s": "s2", "gen_anc": "AFR", "sex": "XY"},
                {"s": "s3", "gen_anc": "EUR", "sex": "XX"},
            ],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

    def test_raw_group_added_by_default(self, strata_table):
        """Test that a 'raw' group (all samples) is inserted at index 1 by default."""
        result = generate_freq_group_membership_array(
            strata_table, [{"gen_anc": strata_table.gen_anc}]
        )
        freq_meta = hl.eval(result.freq_meta)
        counts = hl.eval(result.freq_meta_sample_count)

        # Index 0 is the all-sample adj group, index 1 is the raw group; both
        # cover every sample.
        assert freq_meta[0] == {"group": "adj"}
        assert freq_meta[1] == {"group": "raw"}
        assert counts[0] == counts[1] == strata_table.count()

    def test_no_raw_group_omits_raw(self, strata_table):
        """Test that `no_raw_group=True` omits the raw group entirely."""
        result = generate_freq_group_membership_array(
            strata_table, [{"gen_anc": strata_table.gen_anc}], no_raw_group=True
        )
        freq_meta = hl.eval(result.freq_meta)

        assert freq_meta[0] == {"group": "adj"}
        assert not any(m.get("group") == "raw" for m in freq_meta)

    def test_remove_zero_sample_groups(self, strata_table):
        """Test that `remove_zero_sample_groups` drops groups with no samples."""
        strata = [
            {"gen_anc": strata_table.gen_anc},
            {"sex": strata_table.sex},
            {"gen_anc": strata_table.gen_anc, "sex": strata_table.sex},
        ]

        kept = generate_freq_group_membership_array(strata_table, strata)
        removed = generate_freq_group_membership_array(
            strata_table, strata, remove_zero_sample_groups=True
        )

        # EUR/XY has no samples: kept by default, dropped when requested.
        assert 0 in hl.eval(kept.freq_meta_sample_count)
        assert 0 not in hl.eval(removed.freq_meta_sample_count)
        assert not any(
            m.get("gen_anc") == "EUR" and m.get("sex") == "XY"
            for m in hl.eval(removed.freq_meta)
        )

    def test_multiple_strata_combinations(self, strata_table):
        """Test the cross-product of strata and their sample counts."""
        strata = [
            {"gen_anc": strata_table.gen_anc},
            {"sex": strata_table.sex},
            {"gen_anc": strata_table.gen_anc, "sex": strata_table.sex},
        ]

        result = generate_freq_group_membership_array(strata_table, strata)
        freq_meta = hl.eval(result.freq_meta)
        counts = hl.eval(result.freq_meta_sample_count)

        # Map each stratification group to its sample count (order-independent).
        by_group = {frozenset(m.items()): c for m, c in zip(freq_meta, counts)}
        assert by_group == {
            frozenset({("group", "adj")}): 3,
            frozenset({("group", "raw")}): 3,
            frozenset({("gen_anc", "AFR"), ("group", "adj")}): 2,
            frozenset({("gen_anc", "EUR"), ("group", "adj")}): 1,
            frozenset({("sex", "XX"), ("group", "adj")}): 2,
            frozenset({("sex", "XY"), ("group", "adj")}): 1,
            frozenset({("gen_anc", "AFR"), ("sex", "XX"), ("group", "adj")}): 1,
            frozenset({("gen_anc", "AFR"), ("sex", "XY"), ("group", "adj")}): 1,
            frozenset({("gen_anc", "EUR"), ("sex", "XX"), ("group", "adj")}): 1,
            frozenset({("gen_anc", "EUR"), ("sex", "XY"), ("group", "adj")}): 0,
        }

    @pytest.mark.parametrize(
        "kwargs, match",
        [
            ({"downsamplings": [2]}, "downsampling expression"),
            ({"ds_gen_anc_counts": {"AFR": 1}}, "genetic ancestry group expression"),
        ],
    )
    def test_validation_errors(self, strata_table, kwargs, match):
        """Test that mismatched downsampling/gen-anc inputs raise a ValueError."""
        # Strata lack a `downsampling`/`gen_anc` expression for the given input.
        strata = [{"sex": strata_table.sex}]
        with pytest.raises(ValueError, match=match):
            generate_freq_group_membership_array(strata_table, strata, **kwargs)

    def test_global_downsamplings_without_gen_anc(self):
        """Test downsampling strata without gen_anc use the global ranking.

        A `downsampling` stratum with no `gen_anc` produces global downsampling
        groups whose sample counts equal the requested sizes (`global_idx` < size).
        """
        ht = hl.Table.parallelize(
            [{"s": f"s{i}"} for i in range(6)], hl.tstruct(s=hl.tstr)
        ).key_by("s")
        ds = annotate_downsamplings(ht, [2, 4])

        result = generate_freq_group_membership_array(
            ds, [{"downsampling": ds.downsampling}], downsamplings=[2, 4]
        )
        freq_meta = hl.eval(result.freq_meta)
        counts = hl.eval(result.freq_meta_sample_count)

        by_group = {frozenset(m.items()): c for m, c in zip(freq_meta, counts)}
        # Global downsampling groups hold exactly the requested number of samples.
        assert by_group[frozenset({("downsampling", "2"), ("group", "adj")})] == 2
        assert by_group[frozenset({("downsampling", "4"), ("group", "adj")})] == 4

    def test_group_label_sets_group_key(self):
        """Test that `group_label` sets the 'group' key on every constructed entry."""
        ht = hl.Table.parallelize(
            [{"s": "a", "gen_anc": "AFR"}, {"s": "b", "gen_anc": "EUR"}],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr),
        ).key_by("s")

        result = generate_freq_group_membership_array(
            ht, [{"gen_anc": ht.gen_anc}], group_label="raw", no_raw_group=True
        )
        freq_meta = hl.eval(result.freq_meta)

        assert all(m["group"] == "raw" for m in freq_meta)

    def test_downsampling_group_sample_counts(self):
        """Test that downsampling groups hold the requested number of samples.

        Per-gen_anc downsampling groups select `gen_anc_idx < size` samples from
        the group, and global downsampling groups select `global_idx < size`
        across all samples, so each group's sample count equals the (capped) size.
        """
        rows = [{"s": f"a{i}", "gen_anc": "AFR"} for i in range(5)] + [
            {"s": f"e{i}", "gen_anc": "EUR"} for i in range(4)
        ]
        ht = hl.Table.parallelize(rows, hl.tstruct(s=hl.tstr, gen_anc=hl.tstr)).key_by(
            "s"
        )
        ds = annotate_downsamplings(ht, [3], gen_anc_expr=ht.gen_anc)
        downsamplings = hl.eval(ds.downsamplings)  # [3, 4, 5]
        ds_gen_anc_counts = hl.eval(ds.ds_gen_anc_counts)  # {AFR: 5, EUR: 4}

        result = generate_freq_group_membership_array(
            ds,
            [{"gen_anc": ds.gen_anc, "downsampling": ds.downsampling}],
            downsamplings=downsamplings,
            ds_gen_anc_counts=ds_gen_anc_counts,
        )
        freq_meta = hl.eval(result.freq_meta)
        counts = hl.eval(result.freq_meta_sample_count)
        by_group = {frozenset(m.items()): c for m, c in zip(freq_meta, counts)}

        def ds_count(**kw):
            return by_group[frozenset({**kw, "group": "adj"}.items())]

        # Per-gen_anc: a size below the group count selects exactly that many; a
        # size equal to the group count selects the whole group.
        assert ds_count(gen_anc="AFR", downsampling="3") == 3
        assert ds_count(gen_anc="AFR", downsampling="5") == 5
        assert ds_count(gen_anc="EUR", downsampling="4") == 4
        # Global downsampling groups span all samples.
        assert ds_count(gen_anc="global", downsampling="3") == 3
        assert ds_count(gen_anc="global", downsampling="5") == 5


class TestGksVaFunctions:
    """Tests for GKS/VA helper functions."""

    def test_add_gks_va_grpmax_faf_frequency(self) -> None:
        """
        Test `add_gks_va` output matches expectations for a minimal input struct.

        This test validates:
        - Top-level frequency fields (AC/AN/AF) are propagated from `freq[0]`.
        - QC fields are derived from `filters`, `lcr`, `ab_hist_alt`, and `monoallelic`.
        - `grpMaxFAF95` and `jointGrpMaxFAF95` store numeric FAF values in `frequency`
          (and use the genetic ancestry group label only for the `groupId` suffix).
        """
        # Use a real-looking variant key with a long allele to exercise ID/groupId
        # generation.
        contig = "chr1"
        position = 10108
        ref = "C"
        alt = "CAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT"
        gnomad_id = f"{contig}-{position}-{ref}-{alt}"

        # Construct a minimal Struct that satisfies the fields accessed by `add_gks_va`.
        # Note: keep `filters` empty to avoid nondeterministic ordering from
        # `list(set)`.
        input_struct = hl.eval(
            hl.struct(
                locus=hl.locus(contig, position, reference_genome="GRCh38"),
                alleles=[ref, alt],
                freq=[
                    hl.struct(
                        AC=2,
                        AF=1.98e-04,
                        AN=10090,
                        homozygote_count=0,
                    )
                ],
                filters=hl.empty_set(hl.tstr),
                lcr=False,
                ab_hist_alt=hl.struct(
                    # 20 bins -> 21 edges, matching the expectations in `add_gks_va`.
                    bin_edges=hl.literal(
                        [i / 20 for i in range(21)], hl.tarray(hl.tfloat64)
                    ),
                    # Put counts in the last two bins so the skewed AB count is
                    # non-zero.
                    bin_freq=hl.literal([0] * 18 + [1, 2], hl.tarray(hl.tint64)),
                    n_smaller=hl.int64(0),
                    n_larger=hl.int64(0),
                ),
                in_autosome_or_par=True,
                monoallelic=False,
                grpMaxFAF95=hl.struct(grpmax=7.57e-05, grpmax_gen_anc="nfe"),
                jointGrpMaxFAF95=hl.struct(grpmax=1.23e-04, grpmax_gen_anc="nfe"),
            )
        )

        result = add_gks_va(input_struct=input_struct, label_version="4.1")

        # Validate stable top-level fields derived from the input and label metadata.
        assert result["type"] == "CohortAlleleFrequencyStudyResult"
        assert result["id"] == f"gnomAD-4.1-{gnomad_id}"
        assert result["name"] == f"Overall Cohort Allele Frequency for {gnomad_id}"
        assert result["sourceDataSet"] == {
            "id": "gnomAD4.1",
            "type": "DataSet",
            "name": "gnomAD v4.1",
            "version": "4.1",
        }
        assert result["focusAllele"] == "#/focusAllele"
        assert result["focusAlleleCount"] == 2
        assert result["locusAlleleCount"] == 10090
        assert result["focusAlleleFrequency"] == pytest.approx(1.98e-04)
        assert result["cohort"] == {"id": "ALL", "type": "StudyGroup", "name": "ALL"}

        # Validate ancillary results derived from the input.
        ancillary = result["ancillaryResults"]
        assert ancillary["homozygotes"] == 0
        assert "hemizygotes" not in ancillary
        assert set(ancillary.keys()) == {
            "homozygotes",
            "grpMaxFAF95",
            "jointGrpMaxFAF95",
        }

        # grpMaxFAF95 should store the FAF value and constructed population label.
        grpmax = ancillary["grpMaxFAF95"]
        assert isinstance(grpmax["frequency"], float)
        assert grpmax["frequency"] == pytest.approx(7.57e-05)
        assert grpmax["confidenceInterval"] == 0.95
        assert grpmax["groupId"] == f"{gnomad_id}.NFE"

        # jointGrpMaxFAF95 should behave identically.
        joint_grpmax = ancillary["jointGrpMaxFAF95"]
        assert isinstance(joint_grpmax["frequency"], float)
        assert joint_grpmax["frequency"] == pytest.approx(1.23e-04)
        assert joint_grpmax["confidenceInterval"] == 0.95
        assert joint_grpmax["groupId"] == f"{gnomad_id}.NFE"

        # Validate QC fields derived from filters/LCR/AB histogram/monoallelic.
        quality = result["qualityMeasures"]
        assert quality["qcFilters"] == []
        assert quality["lowComplexityRegion"] is False
        assert quality["heterozygousSkewedAlleleCount"] == 3
        assert quality["monoallelic"] is False
        assert set(quality.keys()) == {
            "qcFilters",
            "lowComplexityRegion",
            "heterozygousSkewedAlleleCount",
            "monoallelic",
        }

        # No subcohort breakdown is expected unless `gen_anc_groups` was requested.
        assert "subCohortFrequency" not in result


class TestCheckAnnotationMissingness:
    """Test the check_annotation_missingness function."""

    @pytest.fixture
    def table_with_simple_annotation(self):
        """Create a Table with a simple annotation containing some missing values."""
        return hl.Table.parallelize(
            [
                {"idx": 1, "value": 10, "optional": "a"},
                {"idx": 2, "value": 20, "optional": None},
                {"idx": 3, "value": None, "optional": "c"},
                {"idx": 4, "value": 40, "optional": None},
                {"idx": 5, "value": 50, "optional": "e"},
            ],
            hl.tstruct(idx=hl.tint32, value=hl.tint32, optional=hl.tstr),
        ).key_by("idx")

    @pytest.fixture
    def table_with_nested_struct(self):
        """Create a Table with nested struct annotations."""
        return hl.Table.parallelize(
            [
                {
                    "idx": 1,
                    "info": hl.Struct(field1=10, field2="a", nested=hl.Struct(x=1.0)),
                },
                {
                    "idx": 2,
                    "info": hl.Struct(field1=20, field2=None, nested=hl.Struct(x=2.0)),
                },
                {
                    "idx": 3,
                    "info": hl.Struct(
                        field1=None, field2="c", nested=hl.Struct(x=None)
                    ),
                },
                {
                    "idx": 4,
                    "info": hl.Struct(field1=40, field2=None, nested=hl.Struct(x=4.0)),
                },
            ],
            hl.tstruct(
                idx=hl.tint32,
                info=hl.tstruct(
                    field1=hl.tint32,
                    field2=hl.tstr,
                    nested=hl.tstruct(x=hl.tfloat64),
                ),
            ),
        ).key_by("idx")

    @pytest.fixture
    def table_with_array_annotation(self):
        """Create a Table with array annotations."""
        return hl.Table.parallelize(
            [
                {"idx": 1, "values": [1, 2, 3]},
                {"idx": 2, "values": []},
                {"idx": 3, "values": None},
                {"idx": 4, "values": [4, 5]},
                {"idx": 5, "values": [None, None]},
            ],
            hl.tstruct(idx=hl.tint32, values=hl.tarray(hl.tint32)),
        ).key_by("idx")

    @pytest.fixture
    def table_with_completely_missing_field(self):
        """Create a Table with a completely missing field."""
        return hl.Table.parallelize(
            [
                {"idx": 1, "present": 10, "missing": None},
                {"idx": 2, "present": 20, "missing": None},
                {"idx": 3, "present": 30, "missing": None},
            ],
            hl.tstruct(idx=hl.tint32, present=hl.tint32, missing=hl.tint32),
        ).key_by("idx")

    @pytest.fixture
    def table_with_multiple_annotations(self):
        """Create a Table with multiple annotations for testing check all."""
        return hl.Table.parallelize(
            [
                {"idx": 1, "field_a": 10, "field_b": "x", "field_c": 1.0},
                {"idx": 2, "field_a": 20, "field_b": None, "field_c": 2.0},
                {"idx": 3, "field_a": None, "field_b": "z", "field_c": None},
                {"idx": 4, "field_a": 40, "field_b": None, "field_c": 4.0},
            ],
            hl.tstruct(
                idx=hl.tint32,
                field_a=hl.tint32,
                field_b=hl.tstr,
                field_c=hl.tfloat64,
            ),
        ).key_by("idx")

    @pytest.fixture
    def matrix_table_with_annotation(self):
        """Create a MatrixTable with row annotations."""
        mt = hl.Table.parallelize(
            [
                {
                    "locus": hl.locus("chr1", 1000, reference_genome="GRCh38"),
                    "alleles": ["A", "T"],
                    "s": "sample1",
                    "GT": hl.call(0, 1),
                },
                {
                    "locus": hl.locus("chr1", 2000, reference_genome="GRCh38"),
                    "alleles": ["C", "G"],
                    "s": "sample1",
                    "GT": hl.call(0, 0),
                },
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])

        # Add row annotations with some missing values.
        mt = mt.annotate_rows(
            info=hl.struct(
                AC=hl.if_else(mt.locus.position == 1000, 5, hl.missing(hl.tint32)),
                AF=0.1,
            )
        )
        return mt

    def test_simple_annotation_missingness(self, table_with_simple_annotation):
        """Test missingness check on simple scalar annotations."""
        ht, results = check_annotation_missingness(
            table_with_simple_annotation, "value"
        )

        assert "missingness" in results
        assert "value" in results["missingness"]
        # 1 out of 5 values is missing = 20%.
        assert results["missingness"]["value"] == pytest.approx(0.2)

    def test_nested_struct_missingness(self, table_with_nested_struct):
        """Test missingness check recursively handles nested structs."""
        ht, results = check_annotation_missingness(table_with_nested_struct, "info")

        assert "missingness" in results
        # Check that all nested fields are present.
        assert "info.field1" in results["missingness"]
        assert "info.field2" in results["missingness"]
        assert "info.nested.x" in results["missingness"]

        # field1: 1 out of 4 missing = 25%.
        assert results["missingness"]["info.field1"] == pytest.approx(0.25)
        # field2: 2 out of 4 missing = 50%.
        assert results["missingness"]["info.field2"] == pytest.approx(0.5)
        # nested.x: 1 out of 4 missing = 25%.
        assert results["missingness"]["info.nested.x"] == pytest.approx(0.25)

    def test_array_annotation_missingness(self, table_with_array_annotation):
        """Test missingness check handles arrays correctly."""
        ht, results = check_annotation_missingness(
            table_with_array_annotation, "values"
        )

        assert "missingness" in results
        assert "values" in results["missingness"]
        # Empty arrays, None, and all-missing-element arrays are considered missing:
        # 3 out of 5 = 60%.
        assert results["missingness"]["values"] == pytest.approx(0.6)

    def test_high_missingness_threshold(self, table_with_nested_struct):
        """Test that fields exceeding high_missingness_threshold are flagged."""
        # With threshold 0.3, field2 (50% missing) should be flagged.
        ht, results = check_annotation_missingness(
            table_with_nested_struct, "info", high_missingness_threshold=0.3
        )

        assert "high_missingness_fields" in results
        assert "info.field2" in results["high_missingness_fields"]
        # field1 and nested.x are at 25%, should not be flagged.
        assert "info.field1" not in results["high_missingness_fields"]
        assert "info.nested.x" not in results["high_missingness_fields"]

    def test_completely_missing_fields_identified(
        self, table_with_completely_missing_field
    ):
        """Test that completely missing fields are identified."""
        ht, results = check_annotation_missingness(
            table_with_completely_missing_field, "missing"
        )

        assert "completely_missing_fields" in results
        assert "missing" in results["completely_missing_fields"]
        assert results["missingness"]["missing"] == 1.0

    def test_remove_missing_fields(self, table_with_completely_missing_field):
        """Test that completely missing fields can be removed."""
        # Add a struct with both present and missing fields.
        ht = table_with_completely_missing_field.annotate(
            info=hl.struct(
                present=table_with_completely_missing_field.present,
                missing=table_with_completely_missing_field.missing,
            )
        )

        ht_result, results = check_annotation_missingness(
            ht, "info", remove_missing_fields=True
        )

        # The 'missing' field should be removed from the struct.
        assert "present" in ht_result.info.dtype.fields
        assert "missing" not in ht_result.info.dtype.fields
        assert "info.missing" in results["completely_missing_fields"]

    def test_matrix_table_row_annotation(self, matrix_table_with_annotation):
        """Test missingness check works on MatrixTable row annotations."""
        mt, results = check_annotation_missingness(matrix_table_with_annotation, "info")

        assert "missingness" in results
        assert "info.AC" in results["missingness"]
        assert "info.AF" in results["missingness"]

        # AC: 1 out of 2 missing = 50%.
        assert results["missingness"]["info.AC"] == pytest.approx(0.5)
        # AF: 0 out of 2 missing = 0%.
        assert results["missingness"]["info.AF"] == pytest.approx(0.0)

    def test_matrix_table_include_col_annotations(self):
        """Test that include_col_annotations checks column annotations."""
        mt = hl.Table.parallelize(
            [
                {
                    "locus": hl.locus("chr1", 1000, reference_genome="GRCh38"),
                    "alleles": ["A", "T"],
                    "s": "sample1",
                    "GT": hl.call(0, 1),
                },
                {
                    "locus": hl.locus("chr1", 1000, reference_genome="GRCh38"),
                    "alleles": ["A", "T"],
                    "s": "sample2",
                    "GT": hl.call(1, 1),
                },
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])

        mt = mt.annotate_cols(
            pop=hl.if_else(mt.s == "sample1", "EUR", hl.missing(hl.tstr))
        )

        mt, results = check_annotation_missingness(mt, include_col_annotations=True)

        # Column annotation should appear with "col." prefix.
        assert "col.pop" in results["missingness"]
        # 1 out of 2 samples missing = 50%.
        assert results["missingness"]["col.pop"] == pytest.approx(0.5)

    def test_annotation_not_found_raises_error(self, table_with_simple_annotation):
        """Test that ValueError is raised when annotation doesn't exist."""
        with pytest.raises(ValueError, match="not found"):
            check_annotation_missingness(
                table_with_simple_annotation, "nonexistent_annotation"
            )

    def test_no_modifications_when_remove_missing_false(
        self, table_with_completely_missing_field
    ):
        """Test that table is not modified when remove_missing_fields=False."""
        ht_result, results = check_annotation_missingness(
            table_with_completely_missing_field,
            "missing",
            remove_missing_fields=False,
        )

        # The field 'missing' should still exist.
        assert "missing" in ht_result.row.dtype.fields
        assert "missing" in results["completely_missing_fields"]

    def test_check_all_annotations_when_none_passed(
        self, table_with_multiple_annotations
    ):
        """Test that all annotations are checked when annotation=None."""
        ht, results = check_annotation_missingness(table_with_multiple_annotations)

        assert "missingness" in results
        # Should have checked all non-key fields.
        assert "field_a" in results["missingness"]
        assert "field_b" in results["missingness"]
        assert "field_c" in results["missingness"]
        # Key field should not be included.
        assert "idx" not in results["missingness"]

        # Verify missingness values.
        # field_a: 1 out of 4 missing = 25%.
        assert results["missingness"]["field_a"] == pytest.approx(0.25)
        # field_b: 2 out of 4 missing = 50%.
        assert results["missingness"]["field_b"] == pytest.approx(0.5)
        # field_c: 1 out of 4 missing = 25%.
        assert results["missingness"]["field_c"] == pytest.approx(0.25)

    def test_check_all_annotations_matrix_table(self, matrix_table_with_annotation):
        """Test checking all annotations on a MatrixTable."""
        mt, results = check_annotation_missingness(matrix_table_with_annotation)

        assert "missingness" in results
        # Should check the info annotation.
        assert "info.AC" in results["missingness"]
        assert "info.AF" in results["missingness"]

    def test_check_all_annotations_with_remove_missing(
        self, table_with_multiple_annotations
    ):
        """Test remove_missing_fields works when checking all annotations."""
        # Add a completely missing field.
        ht = table_with_multiple_annotations.annotate(all_missing=hl.missing(hl.tint32))

        ht_result, results = check_annotation_missingness(
            ht, remove_missing_fields=True
        )

        # The all_missing field should be dropped.
        assert "all_missing" not in ht_result.row.dtype.fields
        assert "all_missing" in results["completely_missing_fields"]
        # Other fields should still exist.
        assert "field_a" in ht_result.row.dtype.fields
        assert "field_b" in ht_result.row.dtype.fields
        assert "field_c" in ht_result.row.dtype.fields


class TestFindMinimalStrataGroups:
    """Test the find_minimal_strata_groups function."""

    def test_basic_adj_with_downsampling(self):
        """Mixed adj/raw partition with gen_anc, sex, and downsampling strata."""
        freq_meta = [
            {"group": "adj"},
            {"group": "raw"},
            {"group": "adj", "gen_anc": "afr"},
            {"group": "adj", "gen_anc": "nfe"},
            {"group": "adj", "gen_anc": "afr", "sex": "XX"},
            {"group": "adj", "gen_anc": "afr", "sex": "XY"},
            {"group": "adj", "gen_anc": "nfe", "sex": "XX"},
            {"group": "adj", "gen_anc": "nfe", "sex": "XY"},
            {"group": "adj", "downsampling": "1000", "gen_anc": "afr"},
            {"group": "adj", "downsampling": "1000", "gen_anc": "nfe"},
            {
                "group": "adj",
                "downsampling": "1000",
                "gen_anc": "afr",
                "sex": "XX",
            },
            {
                "group": "adj",
                "downsampling": "1000",
                "gen_anc": "afr",
                "sex": "XY",
            },
            {
                "group": "adj",
                "downsampling": "1000",
                "gen_anc": "nfe",
                "sex": "XX",
            },
            {
                "group": "adj",
                "downsampling": "1000",
                "gen_anc": "nfe",
                "sex": "XY",
            },
            # Extra raw entries so the raw partition is non-trivial: the
            # bare {"group": "raw"} parent decomposes into the per-gen_anc
            # raw leaves, exercising the adj/raw partition boundary.
            {"group": "raw", "gen_anc": "afr"},
            {"group": "raw", "gen_anc": "nfe"},
        ]
        # 100 samples total: 40 afr (20 XX, 20 XY) + 60 nfe (30 XX, 30 XY).
        # 13 in the downsampling=1000 cohort: 5 afr (2 XX, 3 XY) + 8 nfe (4 XX, 4 XY).
        sample_count = [100, 100, 40, 60, 20, 20, 30, 30, 5, 8, 2, 3, 4, 4, 40, 60]

        leaves, decomp = find_minimal_strata_groups(freq_meta, sample_count)

        assert leaves == [4, 5, 6, 7, 10, 11, 12, 13, 14, 15]
        assert decomp == {
            0: [4, 5, 6, 7],
            1: [14, 15],
            2: [4, 5],
            3: [6, 7],
            8: [10, 11],
            9: [12, 13],
        }

    def test_all_raw_partition(self):
        """compute_stats_per_ref_site case: every entry has group='raw'."""
        freq_meta = [
            {"group": "raw"},
            {"group": "raw", "gen_anc": "afr"},
            {"group": "raw", "gen_anc": "nfe"},
            {"group": "raw", "sex": "XX"},
            {"group": "raw", "sex": "XY"},
            {"group": "raw", "gen_anc": "afr", "sex": "XX"},
            {"group": "raw", "gen_anc": "afr", "sex": "XY"},
            {"group": "raw", "gen_anc": "nfe", "sex": "XX"},
            {"group": "raw", "gen_anc": "nfe", "sex": "XY"},
        ]
        sample_count = [100, 40, 60, 50, 50, 20, 20, 30, 30]

        leaves, decomp = find_minimal_strata_groups(freq_meta, sample_count)

        assert leaves == [5, 6, 7, 8]
        assert decomp == {
            0: [5, 6, 7, 8],
            1: [5, 6],
            2: [7, 8],
            3: [5, 7],
            4: [6, 8],
        }

    def test_force_leaf_groups_retains_pinned_overall_group(self):
        """A group in force_leaf_groups stays a leaf instead of being decomposed.

        Without force_leaf_groups the bare {"group": "raw"} entry (index
        0) is a parent that decomposes into the gen_anc×sex leaves (see
        test_all_raw_partition). Pinning it via force_leaf_groups keeps
        it directly computed — what compute_coverage needs so its
        coverage_stats/qual_hists targets stay cheap index lookups under
        leaf reduction instead of per-row parent reconstruction.
        """
        freq_meta = [
            {"group": "raw"},
            {"group": "raw", "gen_anc": "afr"},
            {"group": "raw", "gen_anc": "nfe"},
            {"group": "raw", "sex": "XX"},
            {"group": "raw", "sex": "XY"},
            {"group": "raw", "gen_anc": "afr", "sex": "XX"},
            {"group": "raw", "gen_anc": "afr", "sex": "XY"},
            {"group": "raw", "gen_anc": "nfe", "sex": "XX"},
            {"group": "raw", "gen_anc": "nfe", "sex": "XY"},
        ]
        sample_count = [100, 40, 60, 50, 50, 20, 20, 30, 30]

        leaves, decomp = find_minimal_strata_groups(
            freq_meta, sample_count, force_leaf_groups=[{"group": "raw"}]
        )

        # Index 0 is now a leaf and absent from the decomposition; the
        # other parents decompose exactly as before.
        assert leaves == [0, 5, 6, 7, 8]
        assert 0 not in decomp
        assert decomp == {
            1: [5, 6],
            2: [7, 8],
            3: [5, 7],
            4: [6, 8],
        }

    def test_global_downsampling_stays_leaf(self):
        """Downsampling-only entries stay leaves; the all-adj entry can't decompose into them."""
        freq_meta = [
            {"group": "adj"},
            {"group": "adj", "downsampling": "1000", "gen_anc": "global"},
            {"group": "adj", "downsampling": "1000", "gen_anc": "afr"},
            {"group": "adj", "downsampling": "1000", "gen_anc": "nfe"},
        ]
        sample_count = [100, 13, 5, 8]

        leaves, decomp = find_minimal_strata_groups(freq_meta, sample_count)

        # Every entry is a leaf:
        # - {} has no matching leaves with non_summable=={}, so falls back to leaf.
        # - downsampling entries are non-summable strata, so they stay as their own leaves.
        assert leaves == [0, 1, 2, 3]
        assert decomp == {}

    def test_no_summable_strata(self):
        """Trivial case: only adj/raw entries, no strata to reduce."""
        freq_meta = [{"group": "adj"}, {"group": "raw"}]
        leaves, decomp = find_minimal_strata_groups(freq_meta, [100, 100])

        assert leaves == [0, 1]
        assert decomp == {}

    def test_partial_strata_promotes_uncoverable_parents_to_leaves(self):
        """When the only candidate decomposition's sample counts don't sum to the parent's, promote to leaf."""
        freq_meta = [
            {"group": "adj"},
            {"group": "adj", "gen_anc": "afr"},
            {"group": "adj", "gen_anc": "afr", "sex": "XX"},
        ]
        # All-adj=100, afr=40, afr-XX=20. The {gen_anc, sex} family is the only
        # leaf family; entry 2 alone doesn't cover entry 0 (20 ≠ 100) or
        # entry 1 (20 ≠ 40), so both must be promoted to leaves.
        sample_count = [100, 40, 20]

        leaves, decomp = find_minimal_strata_groups(freq_meta, sample_count)

        assert leaves == [0, 1, 2]
        assert decomp == {}

    def test_multi_strata_family_uses_sample_count_validation(self):
        """Mirrors the v4 generate_freq pipeline: gen_anc/sex alongside gatk_version/gen_anc.

        Both families have strata set {gen_anc, sex} and {gatk_version, gen_anc}
        respectively — neither is a subset of the other, so both are leaf
        families. Without sample-count validation, parents like the all-adj
        entry would silently sum across both families and double-count
        samples. With validation, each non-leaf is decomposed against exactly
        one family.
        """
        freq_meta = [
            {"group": "adj"},
            {"group": "adj", "gen_anc": "afr"},
            {"group": "adj", "gen_anc": "nfe"},
            {"group": "adj", "sex": "XX"},
            {"group": "adj", "sex": "XY"},
            {"group": "adj", "gen_anc": "afr", "sex": "XX"},
            {"group": "adj", "gen_anc": "afr", "sex": "XY"},
            {"group": "adj", "gen_anc": "nfe", "sex": "XX"},
            {"group": "adj", "gen_anc": "nfe", "sex": "XY"},
            {"group": "adj", "gatk_version": "v1"},
            {"group": "adj", "gatk_version": "v2"},
            {"group": "adj", "gatk_version": "v1", "gen_anc": "afr"},
            {"group": "adj", "gatk_version": "v1", "gen_anc": "nfe"},
            {"group": "adj", "gatk_version": "v2", "gen_anc": "afr"},
            {"group": "adj", "gatk_version": "v2", "gen_anc": "nfe"},
        ]
        # 100 samples: 40 afr (20 XX, 20 XY) + 60 nfe (30 XX, 30 XY).
        # gatk_version split: v1=30 (15 afr + 15 nfe), v2=70 (25 afr + 45 nfe).
        sample_count = [100, 40, 60, 50, 50, 20, 20, 30, 30, 30, 70, 15, 15, 25, 45]

        leaves, decomp = find_minimal_strata_groups(freq_meta, sample_count)

        # Both maximal strata-sets {gen_anc, sex} and {gatk_version, gen_anc}
        # contribute leaves.
        assert leaves == [5, 6, 7, 8, 11, 12, 13, 14]

        # Sex parents (3, 4) only match the {gen_anc, sex} family. gatk
        # parents (9, 10) only match the {gatk_version, gen_anc} family. For
        # parents 0/1/2 both families would sum correctly; the smaller
        # candidate wins, with insertion order breaking ties — the
        # {gen_anc, sex} family is encountered first.
        assert decomp == {
            0: [5, 6, 7, 8],
            1: [5, 6],
            2: [7, 8],
            3: [5, 7],
            4: [6, 8],
            9: [11, 12],
            10: [13, 14],
        }

    def test_fewest_leaves_wins_when_multiple_decompositions_valid(self):
        """When multiple strata families validly cover a parent, the one with the fewest leaves is chosen.

        Two non-comparable leaf families both sum to the parent's sample count:
        `{sex, gatk_version}` with 4 leaves and `{gen_anc}` with 2 leaves. The
        4-leaf family is inserted first, so an insertion-order policy would
        pick it; `min(valid, key=len)` picks the 2-leaf family instead.
        """
        freq_meta = [
            {"group": "adj"},
            # {sex, gatk_version} family — inserted first into the candidate
            # map because its entries have lower indices than the {gen_anc}
            # entries below.
            {"group": "adj", "sex": "XX", "gatk_version": "v1"},
            {"group": "adj", "sex": "XX", "gatk_version": "v2"},
            {"group": "adj", "sex": "XY", "gatk_version": "v1"},
            {"group": "adj", "sex": "XY", "gatk_version": "v2"},
            # {gen_anc} family — fewer leaves, should win.
            {"group": "adj", "gen_anc": "afr"},
            {"group": "adj", "gen_anc": "nfe"},
        ]
        # 100 total. {sex, gatk_version}: 23+27+25+25=100. {gen_anc}: 40+60=100.
        sample_count = [100, 23, 27, 25, 25, 40, 60]

        leaves, decomp = find_minimal_strata_groups(freq_meta, sample_count)

        # Both families are leaves (neither set is a subset of the other).
        assert leaves == [1, 2, 3, 4, 5, 6]
        # The 2-leaf {gen_anc} family wins despite being inserted second.
        assert decomp == {0: [5, 6]}

    def test_custom_non_summable_strata(self):
        """Custom non_summable_strata treats the named stratum as non-summable."""
        freq_meta = [
            {"group": "adj"},
            {"group": "adj", "cohort": "A"},
            {"group": "adj", "cohort": "B"},
        ]
        # 100 total: cohort A=40, cohort B=60.
        sample_count = [100, 40, 60]

        # With cohort treated as a normal summable stratum, entry 0 decomposes
        # into entries 1 and 2 (40+60 == 100).
        leaves, decomp = find_minimal_strata_groups(
            freq_meta, sample_count, non_summable_strata=set()
        )
        assert leaves == [1, 2]
        assert decomp == {0: [1, 2]}

        # With cohort treated as non-summable, no decomposition is possible
        # (entries 1 and 2 are leaves because their non-summable strata don't
        # match entry 0's empty non-summable strata).
        leaves, decomp = find_minimal_strata_groups(
            freq_meta, sample_count, non_summable_strata={"cohort"}
        )
        assert leaves == [0, 1, 2]
        assert decomp == {}

    def test_misaligned_sample_count_raises(self):
        """Length-mismatched sample_count is a programmer error."""
        with pytest.raises(ValueError, match="aligned"):
            find_minimal_strata_groups([{"group": "adj"}], [1, 2])


class TestExpandStrataArrayFromLeaves:
    """Test the expand_strata_array_from_leaves function.

    These cover the per-element-type paths and edge cases that the
    integration tests in `TestAnnotateFreqReduceToMinimalGroups` and
    `TestComputeAlleleNumberPerRefSiteReduceToMinimalGroups` exercise
    only transitively.
    """

    def test_scalar_array_expands_leaves_and_sums_parents(self):
        """Scalar leaf values pass through; parents sum their leaves element-wise."""
        # Original freq_meta has 6 positions; indices 1, 2, 4, 5 are leaves
        # and 0 and 3 are non-leaves with explicit decompositions.
        leaf_array = hl.literal([10, 20, 30, 40], dtype=hl.tarray(hl.tint32))
        leaf_indices = [1, 2, 4, 5]
        decomposition = {0: [1, 2, 4, 5], 3: [4, 5]}
        n_full = 6

        result = hl.eval(
            expand_strata_array_from_leaves(
                leaf_array, leaf_indices, decomposition, n_full
            )
        )
        # 0: all four leaves summed; 1,2,4,5: leaves pass through;
        # 3: leaves at original indices 4 and 5 summed.
        assert result == [100, 10, 20, 70, 30, 40]

    def test_plain_struct_sums_each_field_independently(self):
        """Multi-field structs sum each numeric field independently."""
        element_type = hl.tstruct(x=hl.tint32, y=hl.tint32)
        leaf_array = hl.literal(
            [{"x": 1, "y": 10}, {"x": 2, "y": 20}, {"x": 3, "y": 30}],
            dtype=hl.tarray(element_type),
        )
        leaf_indices = [0, 1, 2]
        decomposition = {3: [0, 1, 2]}
        n_full = 4

        result = hl.eval(
            expand_strata_array_from_leaves(
                leaf_array, leaf_indices, decomposition, n_full
            )
        )
        assert result == [
            hl.Struct(x=1, y=10),
            hl.Struct(x=2, y=20),
            hl.Struct(x=3, y=30),
            hl.Struct(x=6, y=60),
        ]

    def test_freq_struct_drops_and_recomputes_af(self):
        """Freq structs (with AF) drop AF before summing and recompute AF=AC/AN."""
        freq_type = hl.tstruct(
            AC=hl.tint32, AF=hl.tfloat64, AN=hl.tint32, homozygote_count=hl.tint32
        )
        leaf_array = hl.literal(
            [
                {"AC": 2, "AF": 0.04, "AN": 50, "homozygote_count": 0},
                {"AC": 3, "AF": 0.06, "AN": 50, "homozygote_count": 1},
                # Zero-sample leaves: AN=0 must produce missing AF after expand.
                {"AC": 0, "AF": None, "AN": 0, "homozygote_count": 0},
                {"AC": 0, "AF": None, "AN": 0, "homozygote_count": 0},
            ],
            dtype=hl.tarray(freq_type),
        )
        leaf_indices = [1, 2, 4, 5]
        # Parent 0 sums two AN>0 leaves; parent 3 sums two AN=0 leaves.
        decomposition = {0: [1, 2], 3: [4, 5]}
        n_full = 6

        result = hl.eval(
            expand_strata_array_from_leaves(
                leaf_array, leaf_indices, decomposition, n_full
            )
        )
        # Parent of AN>0 leaves: AC=5, AN=100, AF recomputed to 5/100.
        assert result[0] == hl.Struct(AC=5, AF=0.05, AN=100, homozygote_count=1)
        # Leaves with AN>0 are passed through with AF recomputed from AC/AN.
        assert result[1] == hl.Struct(AC=2, AF=0.04, AN=50, homozygote_count=0)
        assert result[2] == hl.Struct(AC=3, AF=0.06, AN=50, homozygote_count=1)
        # Parent of AN=0 leaves: AF must be missing, not divide-by-zero.
        assert result[3] == hl.Struct(AC=0, AF=None, AN=0, homozygote_count=0)
        # Zero-sample leaves: AF stays missing after recompute.
        assert result[4] == hl.Struct(AC=0, AF=None, AN=0, homozygote_count=0)
        assert result[5] == hl.Struct(AC=0, AF=None, AN=0, homozygote_count=0)
        # Field order must match the contract: AC, AF, AN, homozygote_count.
        assert list(result[0].keys()) == ["AC", "AF", "AN", "homozygote_count"]

    def test_missing_leaf_values_treated_as_zero(self):
        """`hl.or_else(..., 0)` zeros out missing leaf values before summing."""
        leaf_array = hl.literal([10, None, 30], dtype=hl.tarray(hl.tint32))
        leaf_indices = [0, 1, 2]
        decomposition = {3: [0, 1, 2]}
        n_full = 4

        result = hl.eval(
            expand_strata_array_from_leaves(
                leaf_array, leaf_indices, decomposition, n_full
            )
        )
        # Leaf position 1 is missing — passes through as 0 because the
        # single-element sum path also applies `or_else(..., 0)`. Parent at
        # position 3 sums 10 + 0 + 30 = 40.
        assert result == [10, 0, 30, 40]

    def test_unaccounted_index_raises(self):
        """Any full-array index that's neither a leaf nor in `decomposition` raises ValueError."""
        leaf_array = hl.literal([1, 2], dtype=hl.tarray(hl.tint32))
        leaf_indices = [0, 1]
        decomposition = {2: [0, 1]}
        # Index 3 is neither in `leaf_indices` nor in `decomposition`.
        n_full = 4

        with pytest.raises(ValueError, match="neither a leaf nor in the decomposition"):
            expand_strata_array_from_leaves(
                leaf_array, leaf_indices, decomposition, n_full
            )


class TestAnnotateFreq:
    """Test the core behavior of annotate_freq.

    Covers the frequency values it produces, the stratifications it builds
    in `freq_meta`/`freq_meta_sample_count`, the `gen_ancs_to_downsample`
    parameter, and the input validation specific to `annotate_freq`.
    """

    @pytest.fixture
    def sample_mt(self):
        """4 samples x 2 variants with hand-computable AC/AN/hom values.

        Samples span two genetic ancestry groups (afr, nfe) and both sexes.
        One genotype (s4 at the first variant) is flagged non-adj so the
        ``adj`` and ``raw`` groups differ.
        """
        samples = [
            ("s1", "afr", "XX"),
            ("s2", "afr", "XY"),
            ("s3", "nfe", "XX"),
            ("s4", "nfe", "XY"),
        ]
        variants = [
            (hl.locus("chr1", 1000, reference_genome="GRCh38"), ["A", "T"]),
            (hl.locus("chr1", 2000, reference_genome="GRCh38"), ["C", "G"]),
        ]
        # Per-variant, per-sample (genotype, adj) in the sample order above.
        cells = [
            # Variant 1: s4 is non-adj so adj != raw here.
            [((0, 1), True), ((1, 1), True), ((0, 0), True), ((0, 1), False)],
            # Variant 2: every genotype is adj.
            [((0, 0), True), ((0, 1), True), ((1, 1), True), ((0, 1), True)],
        ]

        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g, "sex": x} for s, g, x in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

        entries = []
        for v_idx, (locus, alleles) in enumerate(variants):
            for s_idx, (sample_id, _, _) in enumerate(samples):
                (a, b), adj = cells[v_idx][s_idx]
                entries.append(
                    {
                        "locus": locus,
                        "alleles": alleles,
                        "s": sample_id,
                        "GT": hl.call(a, b),
                        "adj": adj,
                    }
                )
        mt = hl.Table.parallelize(
            entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
                adj=hl.tbool,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        return mt.annotate_cols(
            gen_anc=sample_table[mt.s].gen_anc,
            sex=sample_table[mt.s].sex,
        )

    @pytest.fixture
    def gen_anc_downsample_mt(self):
        """7 samples (afr=4, nfe=3) x 1 variant, all adj, with gen_anc annotated.

        The unequal group sizes let downsampling tests distinguish a
        per-group downsampling stratum (afr) from one that is excluded (nfe).
        """
        samples = [("a%i" % i, "afr") for i in range(4)] + [
            ("n%i" % i, "nfe") for i in range(3)
        ]
        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g} for s, g in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr),
        ).key_by("s")
        variant = (hl.locus("chr1", 1000, reference_genome="GRCh38"), ["A", "T"])
        entries = [
            {
                "locus": variant[0],
                "alleles": variant[1],
                "s": s,
                "GT": hl.call(0, 1),
                "adj": True,
            }
            for s, _ in samples
        ]
        mt = hl.Table.parallelize(
            entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
                adj=hl.tbool,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        return mt.annotate_cols(gen_anc=sample_table[mt.s].gen_anc)

    @staticmethod
    def _meta_index(freq_obj):
        """Map each freq_meta entry (as a frozenset of items) to its index."""
        freq_meta = hl.eval(freq_obj.freq_meta)
        return {frozenset(m.items()): i for i, m in enumerate(freq_meta)}

    @staticmethod
    def _freq_by_row(freq_obj):
        """Return the per-variant `freq` arrays keyed by locus position."""
        rows = freq_obj.rows() if isinstance(freq_obj, hl.MatrixTable) else freq_obj
        return {r.locus.position: r.freq for r in rows.select("freq").collect()}

    def test_freq_meta_and_sample_counts(self, sample_mt):
        """freq_meta holds the expected strata with the correct sample counts."""
        mt = annotate_freq(
            sample_mt, sex_expr=sample_mt.sex, gen_anc_expr=sample_mt.gen_anc
        )
        freq_meta = hl.eval(mt.freq_meta)
        counts = hl.eval(mt.freq_meta_sample_count)
        by_group = {frozenset(m.items()): c for m, c in zip(freq_meta, counts)}

        # adj is first, raw is second; both span every sample.
        assert freq_meta[0] == {"group": "adj"}
        assert freq_meta[1] == {"group": "raw"}
        assert by_group == {
            frozenset({("group", "adj")}): 4,
            frozenset({("group", "raw")}): 4,
            frozenset({("group", "adj"), ("sex", "XX")}): 2,
            frozenset({("group", "adj"), ("sex", "XY")}): 2,
            frozenset({("group", "adj"), ("gen_anc", "afr")}): 2,
            frozenset({("group", "adj"), ("gen_anc", "nfe")}): 2,
            frozenset({("group", "adj"), ("gen_anc", "afr"), ("sex", "XX")}): 1,
            frozenset({("group", "adj"), ("gen_anc", "afr"), ("sex", "XY")}): 1,
            frozenset({("group", "adj"), ("gen_anc", "nfe"), ("sex", "XX")}): 1,
            frozenset({("group", "adj"), ("gen_anc", "nfe"), ("sex", "XY")}): 1,
        }

    def test_adj_and_raw_allele_counts(self, sample_mt):
        """The adj and raw call stats match hand-computed values.

        At variant 1, s4's genotype (0/1) is non-adj, so the adj group drops
        that alt allele while the raw group keeps it.
        """
        mt = annotate_freq(
            sample_mt, sex_expr=sample_mt.sex, gen_anc_expr=sample_mt.gen_anc
        )
        idx = self._meta_index(mt)
        freq = self._freq_by_row(mt)
        adj = idx[frozenset({("group", "adj")})]
        raw = idx[frozenset({("group", "raw")})]

        # Variant 1: adj excludes s4's 0/1; raw includes it.
        v1 = freq[1000]
        assert (v1[adj].AC, v1[adj].AN, v1[adj].homozygote_count) == (3, 6, 1)
        assert v1[adj].AF == pytest.approx(3 / 6)
        assert (v1[raw].AC, v1[raw].AN, v1[raw].homozygote_count) == (4, 8, 1)
        assert v1[raw].AF == pytest.approx(4 / 8)

        # Variant 2: all genotypes adj, so adj == raw.
        v2 = freq[2000]
        assert (v2[adj].AC, v2[adj].AN, v2[adj].homozygote_count) == (4, 8, 1)
        assert (v2[raw].AC, v2[raw].AN, v2[raw].homozygote_count) == (4, 8, 1)

    def test_gen_anc_and_sex_stratified_frequencies(self, sample_mt):
        """Per-gen_anc and per-sex call stats match hand-computed values."""
        mt = annotate_freq(
            sample_mt, sex_expr=sample_mt.sex, gen_anc_expr=sample_mt.gen_anc
        )
        idx = self._meta_index(mt)
        freq = self._freq_by_row(mt)

        def stats(pos, **meta):
            s = freq[pos][idx[frozenset({"group": "adj", **meta}.items())]]
            return (s.AC, s.AN, s.homozygote_count)

        # gen_anc strata (adj genotypes only).
        assert stats(1000, gen_anc="afr") == (3, 4, 1)
        assert stats(1000, gen_anc="nfe") == (0, 2, 0)  # s4 non-adj at variant 1
        assert stats(2000, gen_anc="afr") == (1, 4, 0)
        assert stats(2000, gen_anc="nfe") == (3, 4, 1)

        # sex strata (adj genotypes only).
        assert stats(1000, sex="XX") == (1, 4, 0)
        assert stats(1000, sex="XY") == (2, 2, 1)  # only s2; s4 non-adj
        assert stats(2000, sex="XX") == (2, 4, 1)
        assert stats(2000, sex="XY") == (2, 4, 0)

    def test_annotate_mt_false_returns_table(self, sample_mt):
        """`annotate_mt=False` returns a Table carrying the freq annotations."""
        result = annotate_freq(
            sample_mt,
            sex_expr=sample_mt.sex,
            gen_anc_expr=sample_mt.gen_anc,
            annotate_mt=False,
        )
        assert isinstance(result, hl.Table)
        assert "freq" in result.row
        assert "freq_meta" in result.globals
        # The Table-only path produces the same adj call stats as the MT path.
        adj = self._meta_index(result)[frozenset({("group", "adj")})]
        v1 = self._freq_by_row(result)[1000]
        assert (v1[adj].AC, v1[adj].AN) == (3, 6)

    def test_gen_ancs_to_downsample_restricts_downsampling_groups(
        self, gen_anc_downsample_mt
    ):
        """`gen_ancs_to_downsample` limits which per-group downsampling strata appear.

        With afr=4 and nfe=3 samples and a requested downsampling of 2,
        restricting to ['afr'] should produce per-group downsampling strata
        for afr only (sizes 2 and 4); nfe gets no per-group downsampling.
        """
        mt = annotate_freq(
            gen_anc_downsample_mt,
            gen_anc_expr=gen_anc_downsample_mt.gen_anc,
            downsamplings=[2],
            gen_ancs_to_downsample=["afr"],
        )
        freq_meta = hl.eval(mt.freq_meta)

        # Per-group (non-global) downsampling strata that were created.
        per_group = sorted(
            (m["gen_anc"], int(m["downsampling"]))
            for m in freq_meta
            if "downsampling" in m and m.get("gen_anc") not in (None, "global")
        )
        # afr is downsampled to 2 and its full group size 4; nfe is excluded.
        assert per_group == [("afr", 2), ("afr", 4)]
        # nfe still appears as a regular gen_anc stratum, but gets no per-group
        # downsampling stratum since it was excluded from gen_ancs_to_downsample.
        assert not any(
            m.get("gen_anc") == "nfe" and "downsampling" in m for m in freq_meta
        )

    def test_empty_gen_ancs_to_downsample_skips_per_group_downsampling(
        self, gen_anc_downsample_mt
    ):
        """An empty `gen_ancs_to_downsample` yields only global downsamplings.

        Passing `[]` (as opposed to None) means "no per-genetic-ancestry-group
        downsamplings". The function must not fail when no per-group counts are
        produced; it should still generate the global downsampling strata.
        """
        mt = annotate_freq(
            gen_anc_downsample_mt,
            gen_anc_expr=gen_anc_downsample_mt.gen_anc,
            downsamplings=[2],
            gen_ancs_to_downsample=[],
        )
        freq_meta = hl.eval(mt.freq_meta)

        # A global downsampling stratum exists...
        assert any(
            m.get("downsampling") == "2" and m.get("gen_anc") == "global"
            for m in freq_meta
        )
        # ...but no per-genetic-ancestry-group downsampling strata.
        assert not any(
            "downsampling" in m and m.get("gen_anc") not in (None, "global")
            for m in freq_meta
        )

    def test_global_only_downsampling_without_gen_anc(self):
        """`downsamplings` without `gen_anc_expr` produces only global downsamplings.

        This path generates no `ds_gen_anc_counts`, so it exercises the same
        missing-global handling as an empty `gen_ancs_to_downsample`.
        """
        variant = (hl.locus("chr1", 1000, reference_genome="GRCh38"), ["A", "T"])
        entries = [
            {
                "locus": variant[0],
                "alleles": variant[1],
                "s": "s%i" % i,
                "GT": hl.call(0, 1),
                "adj": True,
            }
            for i in range(6)
        ]
        mt = hl.Table.parallelize(
            entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
                adj=hl.tbool,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])

        mt = annotate_freq(mt, downsamplings=[2, 4])
        freq_meta = hl.eval(mt.freq_meta)

        # Global downsampling strata are present; none are per-gen_anc.
        ds_sizes = sorted(
            int(m["downsampling"]) for m in freq_meta if "downsampling" in m
        )
        assert ds_sizes == [2, 4]
        assert all(
            m.get("gen_anc") in (None, "global")
            for m in freq_meta
            if "downsampling" in m
        )

    @pytest.mark.parametrize(
        "kwargs, match",
        [
            # downsampling_expr requires downsamplings.
            (
                {"_use_downsampling_expr": True, "downsamplings": None},
                "requires `downsamplings`",
            ),
            # gen_ancs_to_downsample requires gen_anc_expr.
            (
                {
                    "gen_ancs_to_downsample": ["afr"],
                    "downsamplings": [2],
                    "_use_gen_anc": False,
                },
                "requires `gen_anc_expr`",
            ),
            # gen_ancs_to_downsample is incompatible with a supplied downsampling_expr.
            (
                {
                    "gen_ancs_to_downsample": ["afr"],
                    "downsamplings": [2],
                    "_use_downsampling_expr": True,
                },
                "only used when",
            ),
        ],
    )
    def test_validation_errors(self, sample_mt, kwargs, match):
        """annotate_freq raises ValueError on incompatible argument combinations."""
        kwargs = dict(kwargs)
        use_gen_anc = kwargs.pop("_use_gen_anc", True)
        use_downsampling_expr = kwargs.pop("_use_downsampling_expr", False)

        call_kwargs = {}
        if use_gen_anc:
            call_kwargs["gen_anc_expr"] = sample_mt.gen_anc
        if use_downsampling_expr:
            # A minimal downsampling_expr carrying both indices.
            call_kwargs["downsampling_expr"] = hl.struct(
                global_idx=hl.int32(0), gen_anc_idx=hl.int32(0)
            )
            call_kwargs["ds_gen_anc_counts"] = {"afr": 2}

        with pytest.raises(ValueError, match=match):
            annotate_freq(sample_mt, **call_kwargs, **kwargs)


class TestAnnotateFreqReduceToMinimalGroups:
    """Test that annotate_freq with reduce_to_minimal_groups=True matches the full output."""

    @pytest.fixture
    def sample_mt(self):
        """8 samples × 4 sites with gen_anc, sex, and per-sample-per-site adj flags."""
        samples = [
            ("s1", "afr", "XX"),
            ("s2", "afr", "XY"),
            ("s3", "afr", "XX"),
            ("s4", "afr", "XY"),
            ("s5", "nfe", "XX"),
            ("s6", "nfe", "XY"),
            ("s7", "nfe", "XX"),
            ("s8", "nfe", "XY"),
        ]
        variants = [
            (hl.locus("chr1", 1000, reference_genome="GRCh38"), ["A", "T"]),
            (hl.locus("chr1", 2000, reference_genome="GRCh38"), ["C", "G"]),
            (hl.locus("chr1", 3000, reference_genome="GRCh38"), ["G", "A"]),
            (hl.locus("chr1", 4000, reference_genome="GRCh38"), ["T", "C"]),
        ]
        # Cycle through a few genotype patterns so we get nontrivial AC/AN.
        gt_patterns = [
            [(0, 0), (0, 1), (1, 1), (0, 1), (0, 0), (0, 1), (1, 1), (0, 1)],
            [(0, 1), (1, 1), (0, 0), (0, 1), (0, 1), (1, 1), (0, 0), (0, 1)],
            [(0, 0), (0, 0), (0, 1), (1, 1), (0, 0), (0, 0), (0, 1), (1, 1)],
            [(1, 1), (0, 1), (0, 1), (0, 0), (1, 1), (0, 1), (0, 1), (0, 0)],
        ]

        sample_table = hl.Table.parallelize(
            [{"s": s, "gen_anc": g, "sex": x} for s, g, x in samples],
            hl.tstruct(s=hl.tstr, gen_anc=hl.tstr, sex=hl.tstr),
        ).key_by("s")

        entries = []
        for v_idx, (locus, alleles) in enumerate(variants):
            for s_idx, (sample_id, _, _) in enumerate(samples):
                a, b = gt_patterns[v_idx][s_idx]
                entries.append(
                    {
                        "locus": locus,
                        "alleles": alleles,
                        "s": sample_id,
                        "GT": hl.call(a, b),
                        # Make every other genotype "non-adj" to exercise the
                        # adj/raw distinction.
                        "adj": (s_idx + v_idx) % 2 == 0,
                    }
                )
        mt = hl.Table.parallelize(
            entries,
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                s=hl.tstr,
                GT=hl.tcall,
                adj=hl.tbool,
            ),
        ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])
        mt = mt.annotate_cols(
            gen_anc=sample_table[mt.s].gen_anc,
            sex=sample_table[mt.s].sex,
        )
        return mt

    def _run_and_index_freq(self, mt):
        """Run annotate_freq and return a dict mapping freq_meta key tuples to per-row freq dicts."""
        freq_meta = hl.eval(mt.freq_meta)
        freq_meta_sample_count = hl.eval(mt.freq_meta_sample_count)
        rows = mt.rows().select("freq").collect()
        # Build a stable representation: list of (sorted_meta_items_tuple,
        # [per_row_freq_dict])
        per_meta = {}
        for i, m in enumerate(freq_meta):
            key = tuple(sorted(m.items()))
            per_meta[key] = {
                "sample_count": freq_meta_sample_count[i],
                "freqs": [
                    {
                        "AC": r.freq[i].AC,
                        "AN": r.freq[i].AN,
                        "homozygote_count": r.freq[i].homozygote_count,
                    }
                    for r in rows
                ],
            }
        return per_meta

    def test_reduce_matches_full_with_gen_anc_and_sex(self, sample_mt):
        """annotate_freq(reduce=True) and annotate_freq(reduce=False) match exactly."""
        full_mt = annotate_freq(
            sample_mt,
            sex_expr=sample_mt.sex,
            gen_anc_expr=sample_mt.gen_anc,
        )
        reduced_mt = annotate_freq(
            sample_mt,
            sex_expr=sample_mt.sex,
            gen_anc_expr=sample_mt.gen_anc,
            reduce_to_minimal_groups=True,
        )

        full_indexed = self._run_and_index_freq(full_mt)
        reduced_indexed = self._run_and_index_freq(reduced_mt)

        assert set(full_indexed.keys()) == set(reduced_indexed.keys())
        for key in full_indexed:
            assert (
                full_indexed[key]["sample_count"]
                == reduced_indexed[key]["sample_count"]
            ), key
            assert full_indexed[key]["freqs"] == reduced_indexed[key]["freqs"], key

    def test_reduce_matches_full_with_downsampling(self, sample_mt):
        """annotate_freq(reduce=True) matches the full output when downsamplings are used.

        Downsamplings use a randomized rank, so we pre-compute the downsampling
        annotations once via annotate_downsamplings and pass them as
        downsampling_expr to both annotate_freq calls. This guarantees that
        both runs see exactly the same per-sample downsampling assignments.
        """
        ds_mt = annotate_downsamplings(sample_mt, [4], gen_anc_expr=sample_mt.gen_anc)
        downsamplings = hl.eval(ds_mt.downsamplings)
        ds_gen_anc_counts = hl.eval(ds_mt.ds_gen_anc_counts)

        full_mt = annotate_freq(
            ds_mt,
            sex_expr=ds_mt.sex,
            gen_anc_expr=ds_mt.gen_anc,
            downsamplings=downsamplings,
            downsampling_expr=ds_mt.downsampling,
            ds_gen_anc_counts=ds_gen_anc_counts,
        )
        reduced_mt = annotate_freq(
            ds_mt,
            sex_expr=ds_mt.sex,
            gen_anc_expr=ds_mt.gen_anc,
            downsamplings=downsamplings,
            downsampling_expr=ds_mt.downsampling,
            ds_gen_anc_counts=ds_gen_anc_counts,
            reduce_to_minimal_groups=True,
        )

        full_indexed = self._run_and_index_freq(full_mt)
        reduced_indexed = self._run_and_index_freq(reduced_mt)

        # The freq_meta key sets must match exactly. In particular, the
        # downsampling-stratified entries must survive the reduction as
        # leaves and not be folded into other parents.
        assert set(full_indexed.keys()) == set(reduced_indexed.keys())
        # Sanity check: at least one downsampling-stratified group exists.
        ds_groups = [k for k in full_indexed if any(p[0] == "downsampling" for p in k)]
        assert ds_groups, "expected at least one downsampling-stratified group"

        for key in full_indexed:
            assert (
                full_indexed[key]["sample_count"]
                == reduced_indexed[key]["sample_count"]
            ), key
            assert full_indexed[key]["freqs"] == reduced_indexed[key]["freqs"], key

    def test_reduce_matches_full_with_multi_family_strata(self, sample_mt):
        """annotate_freq(reduce=True) matches the full output when freq_meta has multiple non-comparable strata families.

        Mirrors the v4 generate_freq pipeline shape: a `gen_anc/sex` family
        (from `build_freq_stratification_list`) and a `platform/gen_anc`
        family (from `additional_strata_expr`). Both are maximal-by-inclusion
        leaf-strata sets, and both validly decompose parents like the
        all-adj entry. Without sample-count validation in
        `find_minimal_strata_groups`, the parent would silently sum across
        both families and double-count samples; this test locks that fix in.
        """
        # Mix platforms across gen_anc and sex so neither family perfectly
        # correlates with the other. With this assignment, every parent
        # entry has at least one valid single-family decomposition.
        platforms = {
            "s1": "A",
            "s2": "A",
            "s3": "B",
            "s4": "B",
            "s5": "A",
            "s6": "B",
            "s7": "A",
            "s8": "B",
        }
        platform_table = hl.Table.parallelize(
            [{"s": s, "platform": p} for s, p in platforms.items()],
            hl.tstruct(s=hl.tstr, platform=hl.tstr),
        ).key_by("s")
        mt = sample_mt.annotate_cols(platform=platform_table[sample_mt.s].platform)

        additional = [
            {"platform": mt.platform},
            {"platform": mt.platform, "gen_anc": mt.gen_anc},
        ]
        full_mt = annotate_freq(
            mt,
            sex_expr=mt.sex,
            gen_anc_expr=mt.gen_anc,
            additional_strata_expr=additional,
        )
        reduced_mt = annotate_freq(
            mt,
            sex_expr=mt.sex,
            gen_anc_expr=mt.gen_anc,
            additional_strata_expr=additional,
            reduce_to_minimal_groups=True,
        )

        full_indexed = self._run_and_index_freq(full_mt)
        reduced_indexed = self._run_and_index_freq(reduced_mt)

        assert set(full_indexed.keys()) == set(reduced_indexed.keys())
        # Sanity check: both families are represented, otherwise the
        # multi-family scenario isn't actually being exercised.
        platform_groups = [
            k for k in full_indexed if any(p[0] == "platform" for p in k)
        ]
        assert platform_groups, "expected at least one platform-stratified group"

        for key in full_indexed:
            assert (
                full_indexed[key]["sample_count"]
                == reduced_indexed[key]["sample_count"]
            ), key
            assert full_indexed[key]["freqs"] == reduced_indexed[key]["freqs"], key
