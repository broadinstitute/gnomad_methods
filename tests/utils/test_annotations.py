"""Tests for the annotations utility module."""

from typing import Dict, List

import hail as hl
import pytest

from gnomad.utils.annotations import fill_missing_key_combinations, missing_struct_expr


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

        print(result)
        print(expected)

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
