"""Tests for the annotations utility module."""

from typing import Dict, List

import hail as hl
import pytest

from gnomad.utils.annotations import (
    annotate_downsamplings,
    check_annotation_missingness,
    fill_missing_key_combinations,
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
