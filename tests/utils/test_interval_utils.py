"""Tests for the intervals utility module."""

import logging

import hail as hl
import pytest

from gnomad.utils.intervals import explode_intervals_to_loci


class TestExplodeIntervalsToLoci:
    """Test the explode_intervals_to_loci function."""

    @pytest.fixture
    def sample_interval_table(self):
        """Fixture to create a sample Hail Table with intervals."""
        intervals = [
            hl.interval(
                start=hl.locus("chr1", 100, "GRCh38"),
                end=hl.locus("chr1", 105, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                start=hl.locus("chr2", 200, "GRCh38"),
                end=hl.locus("chr2", 203, "GRCh38"),
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
        return hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 105, "GRCh38"),
            includes_start=True,
            includes_end=True,
        )

    def test_explode_table_with_includes_start_and_end(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test exploding a table with intervals that include both start and end positions."""
        result_ht = explode_intervals_to_loci(
            sample_interval_table, interval_field="interval", keep_intervals=False
        )

        result = result_ht.collect()

        expected_loci_1 = [hl.locus("chr1", pos, "GRCh38") for pos in range(100, 106)]
        expected_loci_2 = [hl.locus("chr2", pos, "GRCh38") for pos in range(200, 203)]

        result_loci = [row.locus for row in result]
        assert len(result_loci) == len(expected_loci_1) + len(expected_loci_2)
        assert all(locus in result_loci for locus in expected_loci_1)
        assert all(locus in result_loci for locus in expected_loci_2)

    def test_explode_table_keep_intervals(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that the original interval field is retained when keep_intervals=True."""
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
        assert "locus" in result_ht.key
        assert all(hasattr(row, "locus") for row in result)

    def test_explode_table_without_keep_intervals(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that the interval field is dropped when keep_intervals=False."""
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

    def test_explode_interval_expression(
        self, sample_interval_expr: hl.expr.IntervalExpression
    ) -> None:
        """Test exploding a single IntervalExpression returns the correct loci."""
        # Explode the interval expression.
        result = explode_intervals_to_loci(sample_interval_expr)

        # Evaluate the result.
        loci = hl.eval(result)

        # Expected loci (chr1:100-105, both inclusive).
        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(100, 106)]

        # Verify the result.
        assert loci == expected_loci

    def test_explode_interval_expression_excludes_start(self) -> None:
        """Test that the start position is excluded when includes_start=False."""
        interval = hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 105, "GRCh38"),
            includes_start=False,
            includes_end=True,
        )

        # Explode the interval expression.
        result = explode_intervals_to_loci(interval)

        # Evaluate the result.
        loci = hl.eval(result)

        # Expected loci (chr1:101-105, start excluded, end included).
        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(101, 106)]

        # Verify the result.
        assert loci == expected_loci

    def test_explode_interval_expression_excludes_end(self) -> None:
        """Test that the end position is excluded when includes_end=False."""
        interval = hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 105, "GRCh38"),
            includes_start=True,
            includes_end=False,
        )

        # Explode the interval expression.
        result = explode_intervals_to_loci(interval)

        # Evaluate the result.
        loci = hl.eval(result)

        # Expected loci (chr1:100-104, start included, end excluded).
        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(100, 105)]

        # Verify the result.
        assert loci == expected_loci

    def test_explode_interval_expression_excludes_both(self) -> None:
        """Test that both endpoints are excluded when includes_start and includes_end are False."""
        interval = hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 105, "GRCh38"),
            includes_start=False,
            includes_end=False,
        )

        # Explode the interval expression.
        result = explode_intervals_to_loci(interval)

        # Evaluate the result.
        loci = hl.eval(result)

        # Expected loci (chr1:101-104, both excluded).
        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(101, 105)]

        # Verify the result.
        assert loci == expected_loci

    def test_explode_table_single_position_interval(self) -> None:
        """Test exploding a table with a single-position interval that includes both endpoints."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 100, "GRCh38"),
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
        assert result[0].locus == hl.locus("chr1", 100, "GRCh38")

    def test_explode_table_missing_interval_field_raises_error(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that a ValueError is raised when interval_field is not provided for a Table."""
        with pytest.raises(
            ValueError, match="`interval_field` and `keep_intervals` must be defined"
        ):
            explode_intervals_to_loci(sample_interval_table)

    def test_explode_table_invalid_interval_field_raises_error(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that an error is raised when interval_field is not present in the Table."""
        with pytest.raises(ValueError, match="`interval_field` must be an annotation"):
            explode_intervals_to_loci(
                sample_interval_table,
                interval_field="nonexistent_field",
                keep_intervals=False,
            )

    def test_explode_table_interval_field_wrong_type_raises_error(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that a ValueError is raised when interval_field refers to a non-interval field."""
        with pytest.raises(ValueError, match="must refer to an interval field"):
            explode_intervals_to_loci(
                sample_interval_table,
                interval_field="gene",
                keep_intervals=False,
            )

    def test_explode_table_grch37(self) -> None:
        """Test exploding a table with GRCh37 reference genome intervals."""
        interval = hl.interval(
            start=hl.locus("1", 1000, "GRCh37"),
            end=hl.locus("1", 1003, "GRCh37"),
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
        expected_loci = [hl.locus("1", pos, "GRCh37") for pos in range(1000, 1004)]

        # Get the loci from the result.
        result_loci = [row.locus for row in result]

        # Verify the result.
        assert len(result_loci) == len(expected_loci)
        assert all(locus in result_loci for locus in expected_loci)

    def test_explode_table_preserves_other_fields(self) -> None:
        """Test that non-interval fields are preserved after exploding."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 102, "GRCh38"),
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
        """Test that an AssertionError is raised when input is an unsupported type."""
        with pytest.raises(
            AssertionError,
            match=(
                "Input must be a Table, IntervalExpression, or list of"
                " IntervalExpressions"
            ),
        ):
            explode_intervals_to_loci("invalid_input")

    def test_explode_interval_expression_single_position_excludes_both(self) -> None:
        """Test that a single-position interval with both endpoints excluded returns an empty array."""
        interval = hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 100, "GRCh38"),
            includes_start=False,
            includes_end=False,
        )

        result = explode_intervals_to_loci(interval)
        loci = hl.eval(result)

        assert loci == []

    def test_explode_table_single_position_interval_excludes_both(self) -> None:
        """Test that a table with a single-position interval excluding both endpoints returns 0 rows."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 100, "GRCh38"),
            includes_start=False,
            includes_end=False,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1"}],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )

        assert result_ht.count() == 0

    def test_explode_table_excludes_start(self) -> None:
        """Test that the start position is excluded when exploding a table interval with includes_start=False."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 105, "GRCh38"),
            includes_start=False,
            includes_end=True,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1"}],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )
        result_loci = [row.locus for row in result_ht.collect()]

        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(101, 106)]
        assert len(result_loci) == len(expected_loci)
        assert all(locus in result_loci for locus in expected_loci)

    def test_explode_table_excludes_end(self) -> None:
        """Test that the end position is excluded when exploding a table interval with includes_end=False."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 105, "GRCh38"),
            includes_start=True,
            includes_end=False,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1"}],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )
        result_loci = [row.locus for row in result_ht.collect()]

        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(100, 105)]
        assert len(result_loci) == len(expected_loci)
        assert all(locus in result_loci for locus in expected_loci)

    def test_explode_table_excludes_both(self) -> None:
        """Test that both endpoints are excluded when exploding a table interval with includes_start and includes_end set to False."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 105, "GRCh38"),
            includes_start=False,
            includes_end=False,
        )

        ht = hl.Table.parallelize(
            [{"interval": interval, "gene": "GENE1"}],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        )
        result_loci = [row.locus for row in result_ht.collect()]

        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(101, 105)]
        assert len(result_loci) == len(expected_loci)
        assert all(locus in result_loci for locus in expected_loci)

    def test_explode_list_of_interval_expressions(self) -> None:
        """Test exploding a list of non-overlapping IntervalExpressions returns all expected loci."""
        intervals = [
            hl.interval(
                hl.locus("chr1", 100, "GRCh38"),
                hl.locus("chr1", 102, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                hl.locus("chr2", 200, "GRCh38"),
                hl.locus("chr2", 202, "GRCh38"),
                includes_start=True,
                includes_end=False,
            ),
        ]

        result = explode_intervals_to_loci(intervals)
        loci = hl.eval(result)

        expected_loci = {
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 101, "GRCh38"),
            hl.locus("chr1", 102, "GRCh38"),
            hl.locus("chr2", 200, "GRCh38"),
            hl.locus("chr2", 201, "GRCh38"),
        }
        assert set(loci) == expected_loci

    def test_explode_table_overlapping_intervals_deduplicates(self) -> None:
        """Test that overlapping table intervals produce deduplicated loci when deduplicate=True."""
        # chr1:100-105 and chr1:103-108, both inclusive → without dedup: 12 rows
        # with deduplicate=True (default): 9 distinct loci [100..108]
        intervals = [
            hl.interval(
                start=hl.locus("chr1", 100, "GRCh38"),
                end=hl.locus("chr1", 105, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                start=hl.locus("chr1", 103, "GRCh38"),
                end=hl.locus("chr1", 108, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
        ]

        ht = hl.Table.parallelize(
            [
                {"interval": intervals[0], "gene": "GENE1"},
                {"interval": intervals[1], "gene": "GENE2"},
            ],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False, deduplicate=True
        )
        result_loci = [row.locus for row in result_ht.collect()]

        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(100, 109)]
        assert len(result_loci) == len(expected_loci)
        assert all(locus in result_loci for locus in expected_loci)

    def test_explode_table_overlapping_intervals_no_deduplicate(self) -> None:
        """Test that overlapping table intervals produce duplicate loci when deduplicate=False."""
        # chr1:100-105 and chr1:103-108, both inclusive.
        # Without deduplication: 6 rows from GENE1 + 6 rows from GENE2 = 12 rows total.
        intervals = [
            hl.interval(
                start=hl.locus("chr1", 100, "GRCh38"),
                end=hl.locus("chr1", 105, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                start=hl.locus("chr1", 103, "GRCh38"),
                end=hl.locus("chr1", 108, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
        ]

        ht = hl.Table.parallelize(
            [
                {"interval": intervals[0], "gene": "GENE1"},
                {"interval": intervals[1], "gene": "GENE2"},
            ],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        result_ht = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False, deduplicate=False
        )
        # Overlapping loci (chr1:103-105) appear twice.
        assert result_ht.count() == 12

    def test_explode_table_overlapping_intervals_keep_intervals_warns(
        self, caplog
    ) -> None:
        """Test that a warning is emitted and duplicates remain when keep_intervals=True with overlapping intervals."""
        # With keep_intervals=True and deduplicate=True, deduplication is skipped and
        # a warning is emitted; overlapping positions appear as duplicate rows.

        intervals = [
            hl.interval(
                start=hl.locus("chr1", 100, "GRCh38"),
                end=hl.locus("chr1", 102, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                start=hl.locus("chr1", 101, "GRCh38"),
                end=hl.locus("chr1", 103, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
        ]

        ht = hl.Table.parallelize(
            [
                {"interval": intervals[0], "gene": "GENE1"},
                {"interval": intervals[1], "gene": "GENE2"},
            ],
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
        )

        with caplog.at_level(logging.WARNING):
            result_ht = explode_intervals_to_loci(
                ht, interval_field="interval", keep_intervals=True, deduplicate=True
            )

        assert (
            "`deduplicate=True` has no effect when `keep_intervals=True`" in caplog.text
        )
        # Duplicate loci are present: positions 101 and 102 each appear twice.
        result_loci = [row.locus for row in result_ht.collect()]
        assert (
            len(result_loci)
            == 6
            # 3 from GENE1 (100,101,102) + 3 from GENE2 (101,102,103), positions
            # 101-102 duplicated
        )

    def test_explode_list_overlapping_intervals_deduplicates(self) -> None:
        """Test that overlapping intervals in a list are deduplicated when deduplicate=True."""
        # Two overlapping intervals; deduplicate=True (default) removes duplicate
        # positions.
        intervals = [
            hl.interval(
                hl.locus("chr1", 100, "GRCh38"),
                hl.locus("chr1", 105, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                hl.locus("chr1", 103, "GRCh38"),
                hl.locus("chr1", 108, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
        ]

        result = explode_intervals_to_loci(intervals, deduplicate=True)
        loci = hl.eval(result)

        expected_loci = {hl.locus("chr1", pos, "GRCh38") for pos in range(100, 109)}
        assert set(loci) == expected_loci

    def test_explode_list_overlapping_intervals_no_deduplicate_warns(
        self, caplog
    ) -> None:
        """Test that a warning is emitted and duplicates are present when deduplicate=False with overlapping list intervals."""
        intervals = [
            hl.interval(
                hl.locus("chr1", 100, "GRCh38"),
                hl.locus("chr1", 102, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                hl.locus("chr1", 101, "GRCh38"),
                hl.locus("chr1", 103, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
        ]

        with caplog.at_level(logging.WARNING):
            result = explode_intervals_to_loci(intervals, deduplicate=False)

        assert (
            "Overlapping intervals in the input list may produce duplicate loci"
            in caplog.text
        )
        loci = hl.eval(result)
        # Loci at chr1:101 and chr1:102 appear twice due to overlap.
        assert len(loci) == 6
        assert loci.count(hl.locus("chr1", 101, "GRCh38")) == 2
        assert loci.count(hl.locus("chr1", 102, "GRCh38")) == 2

    def test_explode_empty_list_raises_type_error(self) -> None:
        """Test that passing an empty Python list raises a TypeError.

        Hail cannot infer the element type of an empty Python list, so
        ``hl.array([])`` (called internally) raises a TypeError. This is distinct
        from passing a Hail-typed empty array expression, which would instead
        fail the input-type assertion because ``ArrayExpression`` is not a
        supported input type.
        """
        with pytest.raises(TypeError):
            explode_intervals_to_loci([])

    def test_explode_list_cross_chromosome_same_position_deduplicates_correctly(
        self,
    ) -> None:
        """Test that loci on different chromosomes with the same integer position are not conflated.

        chr1:100-102 and chr2:100-102 share the same position numbers but are
        distinct loci. Locus-based deduplication must preserve all 6 loci; a
        position-only deduplication strategy would incorrectly collapse them to 3.
        """
        intervals = [
            hl.interval(
                hl.locus("chr1", 100, "GRCh38"),
                hl.locus("chr1", 102, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
            hl.interval(
                hl.locus("chr2", 100, "GRCh38"),
                hl.locus("chr2", 102, "GRCh38"),
                includes_start=True,
                includes_end=True,
            ),
        ]

        result = explode_intervals_to_loci(intervals, deduplicate=True)
        loci = hl.eval(result)

        expected_loci = [hl.locus("chr1", pos, "GRCh38") for pos in range(100, 103)] + [
            hl.locus("chr2", pos, "GRCh38") for pos in range(100, 103)
        ]
        assert len(loci) == 6
        assert set(loci) == set(expected_loci)
