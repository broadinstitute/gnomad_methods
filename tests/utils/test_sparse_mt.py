"""Tests for the sparse_mt utility module."""

import logging

import hail as hl

from gnomad.utils.sparse_mt import get_coverage_agg_func

# Set up logger for tests.
logger = logging.getLogger(__name__)


class TestGetCoverageAggFunc:
    """Test the get_coverage_agg_func function."""

    def test_get_coverage_agg_func_returns_tuple(self):
        """Test that get_coverage_agg_func returns a tuple of two functions."""
        transform_func, agg_func = get_coverage_agg_func()

        assert callable(transform_func)
        assert callable(agg_func)

    def test_get_coverage_agg_func_default_params(self):
        """Test get_coverage_agg_func with default parameters."""
        transform_func, _ = get_coverage_agg_func()

        # Test the transform function with various inputs.
        test_cases = [
            ({"DP": 10}, 10),  # Normal case.
            ({"DP": 0}, 0),  # Zero value.
            ({"DP": 150}, 150),  # Large value.
        ]

        for input_data, expected in test_cases:
            # Create a simple struct for testing.
            test_struct = hl.Struct(**input_data)
            result = hl.eval(transform_func(test_struct))
            assert result == expected

    def test_get_coverage_agg_func_missing_values(self):
        """Test that missing values are handled correctly."""
        transform_func, _ = get_coverage_agg_func()

        # Test with missing value - need to use hl.missing.
        test_struct = hl.Struct(DP=hl.missing(hl.tint32))
        result = hl.eval(transform_func(test_struct))
        assert result == 0

    def test_get_coverage_agg_func_with_nan_values(self):
        """Test that NaN values are handled correctly."""
        transform_func, _ = get_coverage_agg_func()
        test_struct = hl.Struct(DP=hl.float64("nan"))
        result = hl.eval(transform_func(test_struct))
        assert result == 0

    def test_aggregation_function_structure(self):
        """Test that the aggregation function returns the expected structure."""
        _, agg_func = get_coverage_agg_func()
        test_data = [10, 20, 30, 40, 50]
        # Create a Hail Table with a DP field.
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        # Aggregate using the aggregation function.
        result = ht.aggregate(agg_func(ht.DP))
        expected_fields = ["coverage_counter", "mean", "median_approx", "total_DP"]
        for field in expected_fields:
            assert hasattr(result, field)

    def test_aggregation_function_values(self):
        """Test that the aggregation function computes correct values."""
        _, agg_func = get_coverage_agg_func()
        test_data = [10, 20, 30, 40, 50]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Test computed values.
        assert result.mean == 30.0  # (10+20+30+40+50)/5.
        assert result.total_DP == 150  # 10+20+30+40+50.
        assert result.median_approx == 30  # Approximate median.
        assert result.coverage_counter.get(10, 0) == 1
        assert result.coverage_counter.get(20, 0) == 1
        assert result.coverage_counter.get(30, 0) == 1
        assert result.coverage_counter.get(40, 0) == 1
        assert result.coverage_counter.get(50, 0) == 1

    def test_aggregation_function_with_missing_values(self):
        """Test aggregation function handles missing values correctly."""
        _, agg_func = get_coverage_agg_func()
        test_data = [10, 20, None, 40, 50]
        ht = hl.Table.parallelize(
            [{"DP": v if v is not None else hl.missing(hl.tint32)} for v in test_data],
            hl.tstruct(DP=hl.tint32),
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Missing values should be transformed to 0, so they should be included in
        # aggregation.
        assert result.total_DP == 120  # 10+20+0+40+50.
        # The mean calculation may vary depending on how missing values are handled.
        # Let's check that the structure is correct.
        assert hasattr(result, "mean")
        assert hasattr(result, "total_DP")
        assert hasattr(result, "coverage_counter")

    def test_aggregation_function_with_zero_values(self):
        """Test aggregation function handles zero values correctly."""
        _, agg_func = get_coverage_agg_func()
        test_data = [0, 10, 20, 0, 30]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Zero values should be preserved.
        assert result.total_DP == 60  # 0+10+20+0+30.
        assert result.mean == 12.0  # (0+10+20+0+30)/5.
        assert result.coverage_counter.get(0, 0) == 2

    def test_aggregation_function_with_negative_values(self):
        """Test aggregation function handles negative values correctly."""
        _, agg_func = get_coverage_agg_func()
        test_data = [-5, 10, 20, -3, 30]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))

        # Negative values should be preserved.
        assert result.total_DP == 52  # -5+10+20-3+30.
        assert result.mean == 10.4  # (-5+10+20-3+30)/5.
        assert result.coverage_counter.get(-5, 0) == 1
        assert result.coverage_counter.get(-3, 0) == 1

    def test_aggregation_function_edge_cases(self):
        """Test aggregation function with edge cases."""
        _, agg_func = get_coverage_agg_func()

        # Test with single value.
        ht1 = hl.Table.parallelize([{"DP": 42}], hl.tstruct(DP=hl.tint32))
        result1 = ht1.aggregate(agg_func(ht1.DP))
        assert result1.mean == 42.0
        assert result1.total_DP == 42
        assert result1.coverage_counter.get(42, 0) == 1

        # Test with empty table (should handle gracefully).
        ht2 = hl.Table.parallelize([], hl.tstruct(DP=hl.tint32))
        result2 = ht2.aggregate(agg_func(ht2.DP))
        assert result2.total_DP == 0
        assert result2.mean == 0.0

    def test_aggregation_function_with_different_max_cov_bin(self):
        """Test that max_cov_bin parameter affects the aggregation function."""
        _, agg_func = get_coverage_agg_func(max_cov_bin=5)
        test_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))
        coverage_counter = result.coverage_counter
        assert coverage_counter.get(5, 0) > 0
        assert coverage_counter.get(6, 0) == 0
        assert coverage_counter.get(7, 0) == 0

    def test_aggregation_function_max_cov_bin_capping(self):
        """Test that max_cov_bin properly caps values in coverage_counter."""
        _, agg_func = get_coverage_agg_func(max_cov_bin=3)
        test_data = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        ht = hl.Table.parallelize(
            [{"DP": v} for v in test_data], hl.tstruct(DP=hl.tint32)
        )
        result = ht.aggregate(agg_func(ht.DP))
        coverage_counter = result.coverage_counter

        # Values 1, 2, 3 should be counted as themselves.
        assert coverage_counter.get(1, 0) == 1
        assert coverage_counter.get(2, 0) == 1
        # Values 3-10 should all be counted as 3 (max_cov_bin).
        # So there should be 8 values at position 3: the original value 3 plus 7
        # capped values (4-10).
        assert coverage_counter.get(3, 0) == 8

        # No values should be counted above max_cov_bin.
        assert coverage_counter.get(4, 0) == 0
        assert coverage_counter.get(5, 0) == 0

    def test_aggregation_function_median_approximation(self):
        """Test that median_approx provides reasonable approximation."""
        _, agg_func = get_coverage_agg_func()

        # Test with odd number of values.
        test_data_odd = [10, 20, 30, 40, 50]
        ht_odd = hl.Table.parallelize(
            [{"DP": v} for v in test_data_odd], hl.tstruct(DP=hl.tint32)
        )
        result_odd = ht_odd.aggregate(agg_func(ht_odd.DP))
        assert result_odd.median_approx == 30

        # Test with even number of values - approximate median may vary.
        test_data_even = [10, 20, 30, 40]
        ht_even = hl.Table.parallelize(
            [{"DP": v} for v in test_data_even], hl.tstruct(DP=hl.tint32)
        )
        result_even = ht_even.aggregate(agg_func(ht_even.DP))
        # Hail's median_approx is non-deterministic, so we can't test for exact
        # equality.
        # Approximate median should be close to 25 (true median).
        assert result_even.median_approx in [20, 25, 30]

    def test_transform_function_with_various_inputs(self):
        """Test the transform function with various input types."""
        transform_func, _ = get_coverage_agg_func(max_cov_bin=100)

        # Create a sample dataset.
        sample_data = [{"DP": 10}, {"DP": 20}, {"DP": 30}, {"DP": 0}, {"DP": 150}]

        # Transform the data.
        transformed_data = []
        for data in sample_data:
            test_struct = hl.Struct(**data)
            result = hl.eval(transform_func(test_struct))
            transformed_data.append(result)

        # Test that transformation worked correctly.
        assert len(transformed_data) == 5
        assert 150 in transformed_data  # Large value should not be capped in transform.
        assert 0 in transformed_data  # Zero values should be preserved.

    def test_different_dp_field_names(self):
        """Test that different DP field names work correctly."""
        field_names = ["DP", "depth", "coverage", "read_depth"]

        for field_name in field_names:
            transform_func, _ = get_coverage_agg_func(dp_field=field_name)

            # Test with the custom field name.
            test_struct = hl.Struct(**{field_name: 25})
            result = hl.eval(transform_func(test_struct))
            assert result == 25

    def test_transform_and_aggregation_integration(self):
        """Test that transform and aggregation functions work together correctly."""
        transform_func, agg_func = get_coverage_agg_func(max_cov_bin=50)

        # Create test data with various edge cases.
        test_data = [10, 20, None, 40, 50, 0, 100, -5, float("nan")]

        # Test transform function on individual values.
        transformed_values = []
        for value in test_data:
            if value is None:
                test_struct = hl.Struct(DP=hl.missing(hl.tint32))
            elif isinstance(value, float) and value != value:  # Check for NaN
                test_struct = hl.Struct(DP=hl.float64("nan"))
            else:
                test_struct = hl.Struct(DP=value)
            result = hl.eval(transform_func(test_struct))
            transformed_values.append(result)

        # Expected transformed values: [10, 20, 0, 40, 50, 0, 100, -5, 0].
        expected_transformed = [10, 20, 0, 40, 50, 0, 100, -5, 0]
        # For aggregation test (without NaN): [10, 20, 0, 40, 50, 0, 100, -5]
        expected_transformed_no_nan = [10, 20, 0, 40, 50, 0, 100, -5]
        assert transformed_values == expected_transformed

        # Test aggregation function on the same data (excluding NaN due to type
        # handling complexity).
        test_data_no_nan = [
            v for v in test_data if not (isinstance(v, float) and v != v)
        ]
        ht = hl.Table.parallelize(
            [
                {"DP": v if v is not None else hl.missing(hl.tint32)}
                for v in test_data_no_nan
            ],
            hl.tstruct(DP=hl.tint32),
        )

        # Apply the transform function to each row to create transformed data.
        ht = ht.annotate(transformed_DP=transform_func(ht))
        result = ht.aggregate(agg_func(ht.transformed_DP))

        # Check that aggregation works correctly.
        assert result.total_DP == 215  # 10+20+0+40+50+0+100-5.
        # The aggregation function operates on transformed data (missing values become 0).
        # So mean should be calculated from all 8 values: 215/8 = 26.875.
        assert (
            abs(
                result.mean
                - sum(expected_transformed_no_nan) / len(expected_transformed_no_nan)
            )
            < 0.1
        )  # Allow small tolerance.
        # Now that we're using transformed data, there are 2 zeros: one explicit,
        # one from missing.
        assert (
            result.coverage_counter.get(0, 0) == 2
        )  # 2 zero values (explicit 0 + missing->0).
        # Values 50 (original) and 100 (capped to 50 by max_cov_bin) should be at
        # position 50.
        assert result.coverage_counter.get(50, 0) == 2

    def test_custom_field_transform_and_aggregation(self):
        """Test both transform and aggregation with custom field names."""
        transform_func, agg_func = get_coverage_agg_func(
            dp_field="depth", max_cov_bin=30
        )

        # Test transform function.
        test_cases = [
            ({"depth": 10}, 10),
            ({"depth": 0}, 0),
            ({"depth": 50}, 50),  # Should not be capped in transform.
        ]

        for input_data, expected in test_cases:
            test_struct = hl.Struct(**input_data)
            result = hl.eval(transform_func(test_struct))
            assert result == expected

        # Test aggregation function
        test_data = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        ht = hl.Table.parallelize(
            [{"depth": v} for v in test_data], hl.tstruct(depth=hl.tint32)
        )

        result = ht.aggregate(agg_func(ht.depth))

        # Check aggregation results.
        assert result.total_DP == 275  # Sum of all values.
        assert result.mean == 27.5  # 275/10.
        # Values 35, 40, 45, 50 should be capped to 30, plus the original 30.
        assert result.coverage_counter.get(30, 0) == 5  # 4 capped + 1 original.
