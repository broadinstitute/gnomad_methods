"""Tests for the annotations utility module."""

from typing import Dict, List

import hail as hl
import pytest

from gnomad.utils.annotations import (
    fill_missing_key_combinations,
    get_copy_state_by_sex,
    merge_array_expressions,
    merge_freq_arrays,
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
        # Create arrays of integers.
        ht = sample_ht.annotate(
            array1=hl.array([10, 20, 30]),
            array2=hl.array([5, 15, 25]),
        )

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        meta2 = [{"group": "A"}, {"group": "B"}, {"group": "D"}]

        # Merge arrays.
        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2], [meta1, meta2], operation="sum"
        )

        # Collect results.
        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C, D with summed values.
        # A: 10+5=15, B: 20+15=35, C: 30+0=30, D: 0+25=25.
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}, {"group": "D"}]
        assert result_meta == expected_meta

        # Check that arrays are properly merged.
        for row in result:
            assert len(row.result_array) == 4

    def test_merge_integer_arrays_diff(self, sample_ht):
        """Test merging integer arrays with diff operation."""
        # Create arrays of integers.
        ht = sample_ht.annotate(
            array1=hl.array([10, 20, 30]),
            array2=hl.array([5, 15, 25]),
        )

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        meta2 = [{"group": "A"}, {"group": "B"}, {"group": "D"}]

        # Merge arrays with diff operation.
        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2], [meta1, meta2], operation="diff"
        )

        # Collect results.
        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        # Expected: only groups from first array (A, B, C) with subtracted values.
        # A: 10-5=5, B: 20-15=5, C: 30-0=30.
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

    def test_merge_struct_arrays_sum(self, sample_ht):
        """Test merging struct arrays with sum operation."""
        # Create arrays of structs.
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

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge arrays.
        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            struct_fields=["AC", "AN", "homozygote_count"],
        )

        # Collect results.
        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C with summed struct fields.
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        # Check that structs are properly merged.
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
        # Create arrays of structs.
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

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge arrays with diff operation.
        result_array, result_meta = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="diff",
            struct_fields=["AC", "AN", "homozygote_count"],
        )

        # Collect results.
        result = ht.select(result_array=result_array).collect()
        result_meta = hl.eval(result_meta)

        # Expected: only groups from first array (A, B) with subtracted values.
        expected_meta = [{"group": "A"}, {"group": "B"}]
        assert result_meta == expected_meta

        # Check that structs are properly merged.
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
        # Create arrays that will result in negative values.
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([15, 25]),  # Will cause negative values in diff operation.
        )

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Get the result array (this won't raise an error yet).
        result_array, _ = merge_array_expressions(
            [ht.array1, ht.array2], [meta1, meta2], operation="diff"
        )

        # The error should be raised when we try to evaluate the result.
        with pytest.raises(Exception):
            ht.select(result_array=result_array).collect()

    def test_merge_with_negatives_set_to_zero(self, sample_ht):
        """Test that negative values are set to zero when set_negatives_to_zero=True."""
        # Create arrays that will result in negative values.
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([15, 25]),  # Will cause negative values in diff operation.
        )

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Should not raise error and set negatives to zero.
        result_array, _ = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="diff",
            set_negatives_to_zero=True,
        )

        # Collect results.
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
        # Create arrays
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
        )

        # Create count arrays
        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        # Create metadata
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge arrays with count arrays
        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        # Collect results
        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        # Check that both arrays are properly merged
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
        # Create arrays of structs
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

        # Create count arrays
        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        # Create metadata
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge struct arrays with count arrays
        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            struct_fields=["AC", "AN", "homozygote_count"],
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        # Collect results
        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        # Check that both struct and count arrays are properly merged
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
        # Create arrays
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
        )

        # Create multiple count arrays
        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
            qual1=hl.array([30, 40]),
            qual2=hl.array([25, 35]),
        )

        # Create metadata
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge arrays with multiple count arrays
        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="sum",
            count_arrays={
                "sample_count": [ht.count1, ht.count2],
                "quality_score": [ht.qual1, ht.qual2],
            },
        )

        # Collect results
        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        # Check that all arrays are properly merged
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
        # Create arrays
        ht = sample_ht.annotate(
            array1=hl.array([10, 20]),
            array2=hl.array([5, 15]),
        )

        # Create count arrays
        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        # Create metadata
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Merge arrays with count arrays using diff operation
        result_array, result_meta, result_counts = merge_array_expressions(
            [ht.array1, ht.array2],
            [meta1, meta2],
            operation="diff",
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        # Collect results
        result = ht.select(
            result_array=result_array, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        # Expected: only groups from first array (A, B)
        expected_meta = [{"group": "A"}, {"group": "B"}]
        assert result_meta == expected_meta

        # Check that both arrays are properly merged with diff operation
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
        # Create frequency arrays.
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

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge frequency arrays.
        result_freq, result_meta = merge_freq_arrays(
            [ht.freq1, ht.freq2], [meta1, meta2], operation="sum"
        )

        # Collect results.
        result = ht.select(result_freq=result_freq).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C with summed values and calculated AF.
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        # Check that frequency structs are properly merged.
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
        # Create frequency arrays.
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

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge frequency arrays with diff operation.
        result_freq, result_meta = merge_freq_arrays(
            [ht.freq1, ht.freq2], [meta1, meta2], operation="diff"
        )

        # Collect results.
        result = ht.select(result_freq=result_freq).collect()
        result_meta = hl.eval(result_meta)

        # Expected: only groups from first array (A, B) with subtracted values.
        expected_meta = [{"group": "A"}, {"group": "B"}]
        assert result_meta == expected_meta

        # Check that frequency structs are properly merged.
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
        # Create frequency arrays.
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

        # Create count arrays.
        ht = ht.annotate(
            count1=hl.array([100, 200]),
            count2=hl.array([50, 150]),
        )

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "C"}]

        # Merge frequency arrays with count arrays.
        result_freq, result_meta, result_counts = merge_freq_arrays(
            [ht.freq1, ht.freq2],
            [meta1, meta2],
            operation="sum",
            count_arrays={"sample_count": [ht.count1, ht.count2]},
        )

        # Collect results.
        result = ht.select(
            result_freq=result_freq, result_counts=result_counts
        ).collect()
        result_meta = hl.eval(result_meta)

        # Expected: union of groups A, B, C.
        expected_meta = [{"group": "A"}, {"group": "B"}, {"group": "C"}]
        assert result_meta == expected_meta

        # Check that both frequency and count arrays are properly merged.
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

        # Create metadata.
        meta1 = [{"group": "A"}, {"group": "B"}]
        meta2 = [{"group": "A"}, {"group": "B"}]

        # Merge with set_negatives_to_zero=True.
        result_freq, _ = merge_freq_arrays(
            [ht.freq1, ht.freq2],
            [meta1, meta2],
            operation="diff",
            set_negatives_to_zero=True,
        )

        # Collect results.
        result = ht.select(result_freq=result_freq).collect()

        # Check that negative values are set to zero.
        for row in result:
            # A: AC=10-15=-5 -> 0, AN=100-50=50, AF=0/50=0.0.
            first_struct = row.result_freq[0]
            assert first_struct.AC == 0
            assert first_struct.AN == 50
            assert first_struct.AF == 0.0

    def test_merge_freq_arrays_validation_errors(self, sample_ht):
        """Test validation errors in merge_freq_arrays."""
        ht = sample_ht.annotate(
            freq1=hl.array([hl.struct(AC=10, AN=100, homozygote_count=2)]),
            freq2=hl.array([hl.struct(AC=5, AN=50, homozygote_count=1)]),
        )
        meta1 = [{"group": "A"}]
        meta2 = [{"group": "A"}]

        # Test insufficient arrays.
        with pytest.raises(
            ValueError, match="Must provide at least two frequency arrays to merge!"
        ):
            merge_freq_arrays([ht.freq1], [meta1])

        # Test mismatched lengths.
        with pytest.raises(
            ValueError, match="Length of farrays and fmeta must be equal!"
        ):
            merge_freq_arrays([ht.freq1, ht.freq2], [meta1])

        # Test invalid operation.
        with pytest.raises(
            ValueError, match="Operation must be either 'sum' or 'diff'!"
        ):
            merge_freq_arrays([ht.freq1, ht.freq2], [meta1, meta2], operation="invalid")

        # Test mismatched count array lengths.
        with pytest.raises(
            ValueError, match="Length of  count_array 'test' and fmeta must be equal!"
        ):
            merge_freq_arrays(
                [ht.freq1, ht.freq2],
                [meta1, meta2],
                count_arrays={"test": [hl.array([1])]},
            )
