"""Tests for the intervals utility module."""

import logging

import hail as hl
import pytest

from gnomad.utils.intervals import explode_intervals_to_loci

# Endpoint inclusion cases for chr1:100-105.
INCLUSION_CASES = pytest.mark.parametrize(
    "includes_start,includes_end,expected_range",
    [
        (True, True, range(100, 106)),
        (False, True, range(101, 106)),
        (True, False, range(100, 105)),
        (False, False, range(101, 105)),
    ],
    ids=["incl_both", "excl_start", "excl_end", "excl_both"],
)


class TestExplodeIntervalsToLoci:
    """Test the explode_intervals_to_loci function."""

    @pytest.fixture
    def make_interval_table(self):
        """Return a factory that builds a GRCh38 Table with one row per interval."""

        def _make(*intervals):
            return hl.Table.parallelize(
                [
                    {"interval": interval, "gene": f"GENE{i}"}
                    for i, interval in enumerate(intervals, start=1)
                ],
                hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh38")), gene=hl.tstr),
            )

        return _make

    @pytest.fixture
    def sample_interval_table(self, make_interval_table):
        """Fixture to create a sample Hail Table with intervals."""
        return make_interval_table(
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
        )

    @pytest.fixture
    def overlapping_interval_table(self, make_interval_table):
        """Fixture with chr1:100-102 (GENE1) and chr1:101-103 (GENE2), overlapping at 101-102."""
        return make_interval_table(
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

    @pytest.fixture
    def single_pos_exclude_both_interval(self):
        """Single-position interval with both endpoints excluded; represents an empty range."""
        return hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr1", 100, "GRCh38"),
            includes_start=False,
            includes_end=False,
        )

    # --- IntervalExpression input ---

    @INCLUSION_CASES
    def test_explode_interval_expression_inclusion(
        self, sample_interval_expr, includes_start, includes_end, expected_range
    ) -> None:
        """Test that start/end inclusivity is correctly handled for an IntervalExpression."""
        interval = hl.interval(
            sample_interval_expr.start,
            sample_interval_expr.end,
            includes_start=includes_start,
            includes_end=includes_end,
        )
        loci = hl.eval(explode_intervals_to_loci(interval))
        assert loci == [hl.Locus("chr1", pos, "GRCh38") for pos in expected_range]

    def test_explode_interval_expression_single_position_excludes_both(
        self, single_pos_exclude_both_interval
    ) -> None:
        """Test that a single-position interval with both endpoints excluded returns an empty array."""
        assert (
            hl.eval(explode_intervals_to_loci(single_pos_exclude_both_interval)) == []
        )

    # --- Table input ---

    def test_explode_table_with_includes_start_and_end(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test exploding a multi-row table with mixed endpoint inclusion across chromosomes."""
        result_loci = [
            row.locus
            for row in explode_intervals_to_loci(
                sample_interval_table, interval_field="interval", keep_intervals=False
            ).collect()
        ]
        expected = [hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 106)] + [
            hl.Locus("chr2", pos, "GRCh38") for pos in range(200, 203)
        ]
        assert len(result_loci) == len(expected)
        assert all(locus in result_loci for locus in expected)

    def test_explode_table_keep_intervals(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that the original interval field is retained when keep_intervals=True."""
        result_ht = explode_intervals_to_loci(
            sample_interval_table, interval_field="interval", keep_intervals=True
        )
        result = result_ht.collect()
        assert "interval" in result_ht.row
        assert all(hasattr(row, "interval") for row in result)
        assert "locus" in result_ht.key
        assert all(hasattr(row, "locus") for row in result)

    def test_explode_table_without_keep_intervals(
        self, sample_interval_table: hl.Table
    ) -> None:
        """Test that the interval field is dropped when keep_intervals=False."""
        result_ht = explode_intervals_to_loci(
            sample_interval_table, interval_field="interval", keep_intervals=False
        )
        result = result_ht.collect()
        assert "interval" not in result_ht.row
        assert all(not hasattr(row, "interval") for row in result)
        assert all(hasattr(row, "gene") for row in result)

    @INCLUSION_CASES
    def test_explode_table_interval_inclusion(
        self,
        sample_interval_expr,
        make_interval_table,
        includes_start,
        includes_end,
        expected_range,
    ) -> None:
        """Test that start/end inclusivity is correctly handled for Table input."""
        interval = hl.interval(
            sample_interval_expr.start,
            sample_interval_expr.end,
            includes_start=includes_start,
            includes_end=includes_end,
        )
        result_loci = [
            row.locus
            for row in explode_intervals_to_loci(
                make_interval_table(interval),
                interval_field="interval",
                keep_intervals=False,
            ).collect()
        ]
        expected = [hl.Locus("chr1", pos, "GRCh38") for pos in expected_range]
        assert len(result_loci) == len(expected)
        assert all(locus in result_loci for locus in expected)

    def test_explode_table_single_position_interval(self, make_interval_table) -> None:
        """Test exploding a table with a single-position interval that includes both endpoints."""
        interval = hl.interval(
            start=hl.locus("chr1", 100, "GRCh38"),
            end=hl.locus("chr1", 100, "GRCh38"),
            includes_start=True,
            includes_end=True,
        )
        result = explode_intervals_to_loci(
            make_interval_table(interval),
            interval_field="interval",
            keep_intervals=False,
        ).collect()
        assert len(result) == 1
        assert result[0].locus == hl.Locus("chr1", 100, "GRCh38")

    def test_explode_table_single_position_interval_excludes_both(
        self, single_pos_exclude_both_interval, make_interval_table
    ) -> None:
        """Test that a table with a single-position interval excluding both endpoints returns 0 rows."""
        result_ht = explode_intervals_to_loci(
            make_interval_table(single_pos_exclude_both_interval),
            interval_field="interval",
            keep_intervals=False,
        )
        assert result_ht.count() == 0

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
        with pytest.raises(ValueError, match="must be an IntervalExpression"):
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
            hl.tstruct(interval=hl.tinterval(hl.tlocus("GRCh37")), gene=hl.tstr),
        )
        result_loci = [
            row.locus
            for row in explode_intervals_to_loci(
                ht, interval_field="interval", keep_intervals=False
            ).collect()
        ]
        expected = [hl.Locus("1", pos, "GRCh37") for pos in range(1000, 1004)]
        assert len(result_loci) == len(expected)
        assert all(locus in result_loci for locus in expected)

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
        result = explode_intervals_to_loci(
            ht, interval_field="interval", keep_intervals=False
        ).collect()
        assert all(row.gene == "GENE1" for row in result)
        assert all(row.score == 42 for row in result)
        assert len(result) == 3  # positions 100, 101, 102

    def test_explode_table_no_deduplicate_keep_intervals(
        self, overlapping_interval_table: hl.Table
    ) -> None:
        """Test that deduplicate=False with keep_intervals=True preserves duplicate loci with the interval column."""
        result_ht = explode_intervals_to_loci(
            overlapping_interval_table,
            interval_field="interval",
            keep_intervals=True,
            deduplicate=False,
        )
        result = result_ht.collect()
        assert "interval" in result_ht.row
        result_loci = [row.locus for row in result]
        # Overlapping positions 101 and 102 each appear once from GENE1 and once
        # from GENE2.
        assert result_loci.count(hl.Locus("chr1", 101, "GRCh38")) == 2
        assert result_loci.count(hl.Locus("chr1", 102, "GRCh38")) == 2

    def test_explode_table_overlapping_intervals_deduplicates(
        self, overlapping_interval_table: hl.Table
    ) -> None:
        """Test that overlapping table intervals produce deduplicated loci when deduplicate=True."""
        result_loci = [
            row.locus
            for row in explode_intervals_to_loci(
                overlapping_interval_table,
                interval_field="interval",
                keep_intervals=False,
                deduplicate=True,
            ).collect()
        ]
        expected = [hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 104)]
        assert len(result_loci) == len(expected)
        assert all(locus in result_loci for locus in expected)

    def test_explode_table_overlapping_intervals_no_deduplicate(
        self, overlapping_interval_table: hl.Table
    ) -> None:
        """Test that overlapping table intervals produce duplicate loci when deduplicate=False."""
        result_ht = explode_intervals_to_loci(
            overlapping_interval_table,
            interval_field="interval",
            keep_intervals=False,
            deduplicate=False,
        )
        # 3 rows from GENE1 + 3 rows from GENE2.
        assert result_ht.count() == 6

    def test_explode_table_overlapping_intervals_keep_intervals_warns(
        self, overlapping_interval_table: hl.Table, caplog
    ) -> None:
        """Test that a warning is emitted and duplicates remain when keep_intervals=True with overlapping intervals."""
        with caplog.at_level(logging.WARNING):
            result_ht = explode_intervals_to_loci(
                overlapping_interval_table,
                interval_field="interval",
                keep_intervals=True,
                deduplicate=True,
            )
        assert (
            "`deduplicate=True` has no effect when `keep_intervals=True`" in caplog.text
        )
        result_loci = [row.locus for row in result_ht.collect()]
        assert (
            len(result_loci)
            == 6
            # 3 from GENE1 (100,101,102) + 3 from GENE2 (101,102,103)
        )

    # --- ArrayExpression input ---

    def test_explode_invalid_input_type_raises_error(self) -> None:
        """Test that an AssertionError is raised when input is an unsupported type."""
        with pytest.raises(
            AssertionError,
            match=(
                "Input must be a Table, IntervalExpression, or ArrayExpression of"
                " IntervalExpressions"
            ),
        ):
            explode_intervals_to_loci("invalid_input")

    def test_explode_plain_list_raises_assertion_error(self) -> None:
        """Test that passing a plain Python list raises an AssertionError."""
        with pytest.raises(AssertionError):
            explode_intervals_to_loci([])

    def test_explode_list_of_interval_expressions(self) -> None:
        """Test exploding an ArrayExpression of non-overlapping IntervalExpressions returns all expected loci."""
        intervals = hl.array(
            [
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
        )
        loci = hl.eval(explode_intervals_to_loci(intervals))
        expected = {
            hl.Locus("chr1", 100, "GRCh38"),
            hl.Locus("chr1", 101, "GRCh38"),
            hl.Locus("chr1", 102, "GRCh38"),
            hl.Locus("chr2", 200, "GRCh38"),
            hl.Locus("chr2", 201, "GRCh38"),
        }
        assert set(loci) == expected

    def test_explode_array_result_is_sorted(self) -> None:
        """Test that the result is sorted by locus when flatten=True, regardless of input order."""
        # Intervals provided in reverse chromosome order to verify output is sorted.
        intervals = hl.array(
            [
                hl.interval(
                    hl.locus("chr2", 200, "GRCh38"),
                    hl.locus("chr2", 202, "GRCh38"),
                    includes_start=True,
                    includes_end=True,
                ),
                hl.interval(
                    hl.locus("chr1", 100, "GRCh38"),
                    hl.locus("chr1", 102, "GRCh38"),
                    includes_start=True,
                    includes_end=True,
                ),
            ]
        )
        loci = hl.eval(explode_intervals_to_loci(intervals))
        expected = [hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 103)] + [
            hl.Locus("chr2", pos, "GRCh38") for pos in range(200, 203)
        ]
        assert loci == expected

    def test_explode_array_flatten_false(self) -> None:
        """Test that flatten=False returns a nested array of per-interval loci arrays."""
        intervals = hl.array(
            [
                hl.interval(
                    hl.locus("chr1", 100, "GRCh38"),
                    hl.locus("chr1", 102, "GRCh38"),
                    includes_start=True,
                    includes_end=True,
                ),
                hl.interval(
                    hl.locus("chr2", 200, "GRCh38"),
                    hl.locus("chr2", 201, "GRCh38"),
                    includes_start=True,
                    includes_end=False,
                ),
            ]
        )
        result = hl.eval(explode_intervals_to_loci(intervals, flatten=False))
        assert len(result) == 2
        assert result[0] == [hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 103)]
        assert result[1] == [hl.Locus("chr2", 200, "GRCh38")]

    def test_explode_list_overlapping_intervals_deduplicates(self) -> None:
        """Test that overlapping intervals in an ArrayExpression are deduplicated when deduplicate=True."""
        intervals = hl.array(
            [
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
        )
        loci = hl.eval(explode_intervals_to_loci(intervals, deduplicate=True))
        assert set(loci) == {hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 109)}

    def test_explode_list_overlapping_intervals_no_deduplicate_warns(
        self, caplog
    ) -> None:
        """Test that a warning is emitted and duplicates are present when deduplicate=False with overlapping array intervals."""
        intervals = hl.array(
            [
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
        )
        with caplog.at_level(logging.WARNING):
            result = explode_intervals_to_loci(intervals, deduplicate=False)
        assert (
            "Overlapping intervals in the input array may produce duplicate loci"
            in caplog.text
        )
        loci = hl.eval(result)
        assert len(loci) == 6
        assert loci.count(hl.Locus("chr1", 101, "GRCh38")) == 2
        assert loci.count(hl.Locus("chr1", 102, "GRCh38")) == 2

    def test_explode_list_cross_chromosome_independent_intervals_deduplicate_correctly(
        self,
    ) -> None:
        """Test that loci on different chromosomes with the same position number are not conflated.

        chr1:100-102 and chr2:100-102 are separate single-chromosome intervals sharing
        the same position numbers. Locus-based deduplication must preserve all 6 loci.
        """
        intervals = hl.array(
            [
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
        )
        loci = hl.eval(explode_intervals_to_loci(intervals, deduplicate=True))
        expected = [hl.Locus("chr1", pos, "GRCh38") for pos in range(100, 103)] + [
            hl.Locus("chr2", pos, "GRCh38") for pos in range(100, 103)
        ]
        assert len(loci) == 6
        assert set(loci) == set(expected)

    def test_explode_cross_chromosome_interval_raises_error(self) -> None:
        """Test that a cross-chromosome interval raises an error at evaluation time."""
        cross_chrom = hl.interval(
            hl.locus("chr1", 100, "GRCh38"),
            hl.locus("chr2", 200, "GRCh38"),
        )
        result = explode_intervals_to_loci(cross_chrom)
        with pytest.raises(Exception, match="[Cc]ross-chromosome"):
            hl.eval(result)
