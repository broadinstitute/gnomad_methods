"""Tests for the file_utils utility module."""

import logging

import hail as hl

from gnomad.utils.file_utils import (
    convert_multi_array_to_array_of_structs,
    print_global_struct,
)


class TestPrintGlobalStruct:
    """Test the print_global_struct function."""

    def test_with_table(self, caplog):
        """Test output format when passing a Table."""
        ht = hl.Table.parallelize(
            [{"x": 1}],
            hl.tstruct(x=hl.tint32),
        )
        ht = ht.annotate_globals(foo="bar", count=42)

        with caplog.at_level(logging.INFO):
            print_global_struct(ht)

        assert "foo: bar" in caplog.text
        assert "count: 42" in caplog.text

    def test_with_struct(self, caplog):
        """Test output format with an already-evaluated Struct."""
        s = hl.Struct(a=1, b="hello", nested=hl.Struct(c=3.0))

        with caplog.at_level(logging.INFO):
            print_global_struct(s)

        assert "a: 1" in caplog.text
        assert "b: hello" in caplog.text
        assert "c: 3.0" in caplog.text

    def test_with_struct_expression(self, caplog):
        """Test output format when passing a StructExpression."""
        ht = hl.Table.parallelize(
            [{"x": 1}],
            hl.tstruct(x=hl.tint32),
        )
        ht = ht.annotate_globals(foo="bar", count=42)

        with caplog.at_level(logging.INFO):
            print_global_struct(ht.globals)

        assert "foo: bar" in caplog.text
        assert "count: 42" in caplog.text

    def test_nested_indentation(self, caplog):
        """Test that nested structs are indented deeper than top-level fields."""
        s = hl.Struct(top="val", nested=hl.Struct(inner="deep"))

        with caplog.at_level(logging.INFO):
            print_global_struct(s)

        # Top-level fields get 4 spaces, nested get 8.
        assert "    top: val" in caplog.text
        assert "        inner: deep" in caplog.text

    def test_multiple_nested_levels(self, caplog):
        """Test formatting with multiple nesting levels."""
        s = hl.Struct(level1=hl.Struct(level2=hl.Struct(value=99)))

        with caplog.at_level(logging.INFO):
            print_global_struct(s)

        assert "    level1:" in caplog.text
        assert "        level2:" in caplog.text
        assert "            value: 99" in caplog.text


class TestConvertMultiArrayToArrayOfStructs:
    """Test the convert_multi_array_to_array_of_structs function."""

    def test_basic_conversion(self):
        """Test converting two parallel arrays into an array of structs."""
        ht = hl.Table.parallelize(
            [{"a": [1, 2, 3], "b": [4, 5, 6], "other": "keep"}],
            hl.tstruct(
                a=hl.tarray(hl.tint32),
                b=hl.tarray(hl.tint32),
                other=hl.tstr,
            ),
        )

        result_ht = convert_multi_array_to_array_of_structs(ht, ["a", "b"], "combined")
        result = result_ht.collect()[0]

        # Original fields should be dropped.
        assert not hasattr(result, "a")
        assert not hasattr(result, "b")

        # Non-combined field should remain.
        assert result.other == "keep"

        # Check the combined array.
        assert len(result.combined) == 3
        assert result.combined[0].a == 1
        assert result.combined[0].b == 4
        assert result.combined[1].a == 2
        assert result.combined[1].b == 5
        assert result.combined[2].a == 3
        assert result.combined[2].b == 6

    def test_three_arrays(self):
        """Test converting three parallel arrays."""
        ht = hl.Table.parallelize(
            [{"x": [1.0, 2.0], "y": [3.0, 4.0], "z": [5.0, 6.0]}],
            hl.tstruct(
                x=hl.tarray(hl.tfloat64),
                y=hl.tarray(hl.tfloat64),
                z=hl.tarray(hl.tfloat64),
            ),
        )

        result_ht = convert_multi_array_to_array_of_structs(
            ht, ["x", "y", "z"], "merged"
        )
        result = result_ht.collect()[0]

        assert len(result.merged) == 2
        assert result.merged[0].x == 1.0
        assert result.merged[0].y == 3.0
        assert result.merged[0].z == 5.0

    def test_with_struct_expression(self):
        """Test conversion on a StructExpression (annotated within a table)."""
        ht = hl.Table.parallelize(
            [{"s": hl.Struct(a=[10, 20], b=[30, 40])}],
            hl.tstruct(
                s=hl.tstruct(
                    a=hl.tarray(hl.tint32),
                    b=hl.tarray(hl.tint32),
                )
            ),
        )

        result_ht = ht.annotate(
            s=convert_multi_array_to_array_of_structs(ht.s, ["a", "b"], "combined")
        )
        result = result_ht.collect()[0]

        assert len(result.s.combined) == 2
        assert result.s.combined[0].a == 10
        assert result.s.combined[0].b == 30
