"""Tests for the annotations utility module."""

from typing import Dict, List

import ga4gh.core as ga4gh_core
import ga4gh.vrs as ga4gh_vrs
import hail as hl
import pytest

from gnomad.utils.annotations import (
    VRS_CHROM_IDS,
    add_gks_vrs,
    annotate_downsamplings,
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
        ], (
            f"Failed for locus={locus}, is_xx={is_xx}. Expected {[expected_diploid, expected_hemi_x, expected_hemi_y]}, got {result}"
        )


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

    def test_vrs_identifier_generation(self):
        """Test that GA4GH identifiers are generated correctly."""
        # Test that we can import and access VRS modules
        assert hasattr(ga4gh_core, "__version__")
        assert hasattr(ga4gh_vrs, "models")
        assert hasattr(ga4gh_vrs.models, "SequenceLocation")

        # Test that VRS 2.0.1+ API is available
        assert hasattr(ga4gh_core, "ga4gh_identify"), (
            "VRS 2.0.1+ ga4gh_identify function not found"
        )

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
        """Test that add_gks_vrs actually calls the VRS 2.0.1 API with real data."""
        # Create a real Hail locus and VRS struct that would come from actual data
        locus = hl.locus("chr1", 100, reference_genome="GRCh38")

        # Create a VRS struct that mimics what would be in actual gnomAD data
        vrs_struct = hl.struct(
            VRS_Allele_IDs=["test_ref_id", "ga4gh:VA.test_var_id"],
            VRS_Starts=[99, 99],  # 0-based coordinates
            VRS_Ends=[100, 100],
            VRS_States=["A", "T"],  # ref, alt
        )

        # Evaluate to get actual Python objects
        locus_py = hl.eval(locus)
        vrs_py = hl.eval(vrs_struct)

        result = add_gks_vrs(locus_py, vrs_py)
        # Validate with Pydantic, but do not renormalize.
        _ = ga4gh_vrs.models.Allele(**result)

        # Verify the result has the expected VRS structure
        assert isinstance(result, dict)
        assert result["type"] == "Allele"
        assert "location" in result
        assert result["location"]["type"] == "SequenceLocation"
        assert "id" in result["location"]  # This comes from ga4gh_identify() call
        assert "state" in result
        assert result["state"]["type"] == "LiteralSequenceExpression"
        assert result["state"]["sequence"] == "T"

        # Verify the chromosome ID matches our VRS_CHROM_IDS mapping
        expected_chr1_id = VRS_CHROM_IDS["GRCh38"]["chr1"]
        assert expected_chr1_id == "SQ.Ya6Rs7DHhDeg7YaOSg1EoNi3U_nQ9SvO"
        assert (
            result["location"]["sequenceReference"]["refgetAccession"]
            == expected_chr1_id
        )

        # Verify the location ID was generated by the VRS API
        location_id = result["location"]["id"]
        assert isinstance(location_id, str)
        assert location_id.startswith("ga4gh:")  # VRS identifiers start with ga4gh:

    def test_add_gks_vrs_grch37_chromosome_mapping(self):
        """Test that VRS chromosome mapping works correctly for GRCh37."""
        # Test GRCh37 (note different chromosome naming: "1" vs "chr1")
        locus_grch37 = hl.locus("1", 100, reference_genome="GRCh37")
        vrs_struct = hl.struct(
            VRS_Allele_IDs=["test_ref_id", "ga4gh:VA.test_var_id"],
            VRS_Starts=[99, 99],
            VRS_Ends=[100, 100],
            VRS_States=["A", "T"],
        )

        locus_py = hl.eval(locus_grch37)
        vrs_py = hl.eval(vrs_struct)

        result = add_gks_vrs(locus_py, vrs_py)
        # Validate with Pydantic, but do not renormalize.
        _ = ga4gh_vrs.models.Allele(**result)

        # Verify GRCh37 chromosome 1 mapping
        expected_chr1_grch37_id = VRS_CHROM_IDS["GRCh37"]["1"]
        assert expected_chr1_grch37_id == "SQ.S_KjnFVz-FE7M0W6yoaUDgYxLPc1jyWU"
        assert (
            result["location"]["sequenceReference"]["refgetAccession"]
            == expected_chr1_grch37_id
        )

        # Verify the API was called and generated a proper identifier
        location_id = result["location"]["id"]
        assert isinstance(location_id, str)
        assert location_id.startswith("ga4gh:")

    def test_vrs_version_compatibility(self):
        """Test that VRS 2.0.1+ is properly installed."""
        # Test that VRS 2.0.1+ API is available
        assert hasattr(ga4gh_core, "ga4gh_identify"), "VRS 2.0.1+ is required"

        # The function should work with VRS 2.0.1+
        assert callable(add_gks_vrs)


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
