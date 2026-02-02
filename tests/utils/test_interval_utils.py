"""Tests for the intervals utility module."""

import hail as hl
import pytest

from gnomad.utils.intervals import explode_intervals_to_loci


class TestExplodeIntervalsToLoci:
    """Test the explode_intervals_to_loci function."""

    @pytest.fixture
    def sample_interval_table(self):
        """Fixture to create a sample Hail Table with intervals."""
        intervals = [
            hl.Interval(
                start=hl.Locus("chr1", 100, "GRCh38"),
                end=hl.Locus("chr1", 105, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.Interval(
                start=hl.Locus("chr2", 200, "GRCh38"),
                end=hl.Locus("chr2", 203, "GRCh38"),
                includes_start=True,
                includes_end=False,
            ),
        ]
        return hl.Table.parallelize(
            [
                {"interval": intervals[0], "gene": "GENE1"},
                {"interval": intervals[1], "gene": "GENE2"},
            ],
            hl.tstruct(
                interval=hl.tinterval(hl.tlocus("GRCh38")),
                gene=hl.tstr,
            ),
        )

    @pytest.fixture
    def sample_interval_expr(self):
        """Fixture to create a sample interval expression."""
        return hl.Interval(
            start=hl.Locus("chr1", 100, "GRCh38"),
            end=hl.Locus("chr1", 105, "GRCh38"),
            includes_start=True,
            includes_end=True,
        )

    def test_explode_table_with_includes_start_and_end(
        self, sample_interval_table: hl.Table
    ) -> None:
        """
        Test exploding a table with intervals that include both start and end positions.

        :param sample_interval_table: Sample Hail Table with intervals.
        :return: None.
        """
        # Explode the intervals to loci.
        result_ht = explode_intervals_to_loci(
            sample_interval_table, interval_field="interval", keep_intervals=False
        )

        # Collect results.
        result = result_ht.collect()

        # Expected loci for the first interval (100-105, both inclusive).
        expected_loci_1 = [
            hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 106)
        ]
        # Expected loci for the second interval (200-203, start inclusive, end exclusive).
        expected_loci_2 = [
            hl.Locus("chr2", pos, "GRCh38") for pos in range(200, 203)
        ]

        # Get the loci from the result.
        result_loci = [row.locus for row in result]

        # Verify the result contains the expected loci.
        assert len(result_loci) == len(expected_loci_1) + len(expected_loci_2)
        assert all(locus in result_loci for locus in expected_loci_1)
        assert all(locus in result_loci for locus in expected_loci_2)

    def test_explode_table_keep_intervals(
        self, sample_interval_table: hl.Table
    ) -> None:
        """
        Test exploding a table while keeping the original interval field.

        :param sample_interval_table: Sample Hail Table with intervals.
        :return: None.
        """
        # Explode the intervals to loci, keeping the interval field.
        result_ht = explode_intervals_to_loci(
            sample_interval_table, interval_field="interval", keep_intervals=True
        )

        # Collect results.
        result = result_ht.collect()

        # Verify that the interval field is still present.
        assert "interval" in result_ht.row
        assert all(hasattr(row, "interval") for row in result)

        # Verify that the locus field is present.
        assert "locus" in result_ht.row_key
        assert all(hasattr(row, "locus") for row in result)

    def test_explode_table_without_keep_intervals(
        self, sample_interval_table: hl.Table
    ) -> None:
        """
        Test exploding a table without keeping the original interval field.

        :param sample_interval_table: Sample Hail Table with intervals.
        :return: None.
        """
        # Explode the intervals to loci without keeping the interval field.
        result_ht = explode_intervals_to_loci(
            sample_interval_table, interval_field="interval", keep_intervals=False
        )

        # Collect results.
        result = result_ht.collect()

        # Verify that the interval field is not present.
        assert "interval" not in result_ht.row
        assert all(not hasattr(row, "interval") for row in result)

        # Verify that other fields are still present.
        assert all(hasattr(row, "gene") for row in result)

    def test_explode_interval_expression(self, sample_interval_expr: hl.Interval) -> None:
        """
        Test exploding an interval expression to an array of positions.

        :param sample_interval_expr: Sample interval expression.
        :return: None.
        """
        # Explode the interval expression.
        result = explode_intervals_to_loci(sample_interval_expr)

        # Evaluate the result.
        positions = hl.eval(result)

        # Expected positions (100-105, both inclusive).
        expected_positions = list(range(100, 106))

        # Verify the result.
        assert positions == expected_positions

    def test_explode_interval_expression_excludes_start(self) -> None:
        """
        Test exploding an interval expression that excludes the start position.

        :return: None.
        """
        interval = hl.Interval(
            start=hl.Locus("chr1", 100, "GRCh38"),
            end=hl.Locus("chr1", 105, "GRCh38"),
            includes_start=False,
            includes_end=True,
        )

        # Explode the interval expression.
        result = explode_intervals_to_loci(interval)

        # Evaluate the result.
        positions = hl.eval(result)

        # Expected positions (101-105, start excluded, end included).
        expected_positions = list(range(101, 106))

        # Verify the result.
        assert positions == expected_positions

    def test_explode_interval_expression_excludes_end(self) -> None:
        """
        Test exploding an interval expression that excludes the end position.

        :return: None.
        """
        interval = hl.Interval(
            start=hl.Locus("chr1", 100, "GRCh38"),
            end=hl.Locus("chr1", 105, "GRCh38"),
            includes_start=True,
            includes_end=False,
        )

        # Explode the interval expression.
        result = explode_intervals_to_loci(interval)

        # Evaluate the result.
        positions = hl.eval(result)

        # Expected positions (100-104, start included, end excluded).
        expected_positions = list(range(100, 105))

        # Verify the result.
        assert positions == expected_positions

    def test_explode_interval_expression_excludes_both(self) -> None:
        """
        Test exploding an interval expression that excludes both start and end positions.

        :return: None.
        """
        interval = hl.Interval(
            start=hl.Locus("chr1", 100, "GRCh38"),
            end=hl.Locus("chr1", 105, "GRCh38"),
            includes_start=False,
            includes_end=False,
        )

        # Explode the interval expression.
        result = explode_intervals_to_loci(interval)

        # Evaluate the result.
        positions = hl.eval(result)

        # Expected positions (101-104, both excluded).
        expected_positions = list(range(101, 105))

        # Verify the result.
        assert positions == expected_positions

    def test_explode_table_single_position_interval(self) -> None:
        """
        Test exploding a table with a single-position interval.

        :return: None.
        """
        interval = hl.Interval(
            start=hl.Locus("chr1", 100, "GRCh38"),
            end=hl.Locus("chr1", 100, "GRCh38"),
            includes_start=True,
            includes_end=True,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1"}],
            hl.tstruct(
                interval=hl.tinterval(hl.tlocus("GRCh38")),
                gene=hl.tstr,
            ),
        )

        # Explode the intervals to loci.
        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )

        # Collect results.
        result = result_ht.collect()

        # Expected: single locus at position 100.
        assert len(result) == 1
        assert result[0].locus == hl.Locus("chr1", 100, "GRCh38")

    def test_explode_table_missing_interval_field_raises_error(
        self, sample_interval_table: hl.Table
    ) -> None:
        """
        Test that an error is raised when interval_field is not provided for a Table.

        :param sample_interval_table: Sample Hail Table with intervals.
        :return: None.
        """
        with pytest.raises(ValueError, match="`interval_field` and `keep_intervals` must be defined"):
            explode_intervals_to_loci(sample_interval_table)

    def test_explode_table_invalid_interval_field_raises_error(
        self, sample_interval_table: hl.Table
    ) -> None:
        """
        Test that an error is raised when interval_field is not in the Table.

        :param sample_interval_table: Sample Hail Table with intervals.
        :return: None.
        """
        with pytest.raises(AssertionError, match="`interval_field` must be an annotation"):
            explode_intervals_to_loci(
                sample_interval_table,
                interval_field="nonexistent_field",
                keep_intervals=False,
            )

    def test_explode_table_grch37(self) -> None:
        """
        Test exploding a table with GRCh37 reference genome.

        :return: None.
        """
        interval = hl.Interval(
            start=hl.Locus("1", 1000, "GRCh37"),
            end=hl.Locus("1", 1003, "GRCh37"),
            includes_start=True,
            includes_end=True,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1"}],
            hl.tstruct(
                interval=hl.tinterval(hl.tlocus("GRCh37")),
                gene=hl.tstr,
            ),
        )

        # Explode the intervals to loci.
        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )

        # Collect results.
        result = result_ht.collect()

        # Expected loci (1000-1003, both inclusive).
        expected_loci = [hl.Locus("1", pos, "GRCh37") for pos in range(1000, 1004)]

        # Get the loci from the result.
        result_loci = [row.locus for row in result]

        # Verify the result.
        assert len(result_loci) == len(expected_loci)
        assert all(locus in result_loci for locus in expected_loci)

    def test_explode_table_preserves_other_fields(self) -> None:
        """
        Test that exploding preserves other fields in the table.

        :return: None.
        """
        interval = hl.Interval(
            start=hl.Locus("chr1", 100, "GRCh38"),
            end=hl.Locus("chr1", 102, "GRCh38"),
            includes_start=True,
            includes_end=True,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1", "score": 42}],
            hl.tstruct(
                interval=hl.tinterval(hl.tlocus("GRCh38")),
                gene=hl.tstr,
                score=hl.tint32,
            ),
        )

        # Explode the intervals to loci.
        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )

        # Collect results.
        result = result_ht.collect()

        # Verify that other fields are preserved.
        assert all(row.gene == "GENE1" for row in result)
        assert all(row.score == 42 for row in result)
        assert len(result) == 3  # 3 positions: 100, 101, 102

    def test_explode_invalid_input_type_raises_error(self) -> None:
        """
        Test that an error is raised when input is neither a Table nor an IntervalExpression.

        :return: None.
        """
        with pytest.raises(AssertionError, match="Input must be a Table or IntervalExpression"):
            explode_intervals_to_loci("invalid_input")
