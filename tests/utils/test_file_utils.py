"""Tests for the file_utils utility module."""

import logging

import hail as hl
import pytest
from hail.genetics import Locus
from hail.utils import Interval, Struct

from gnomad.utils.file_utils import repartition_for_join

# Set up logger for tests.
logger = logging.getLogger(__name__)


@pytest.fixture(scope="module")
def int_keyed_ht():
    """Return a range Table keyed by `idx` with 10 partitions."""
    return hl.utils.range_table(100, n_partitions=10)


@pytest.fixture(scope="module")
def locus_keyed_ht():
    """Return a Table keyed by `locus` with 10 partitions."""
    ht = hl.utils.range_table(100, n_partitions=10)
    ht = ht.key_by(locus=hl.locus("chr1", ht.idx + 1, reference_genome="GRCh38"))
    return ht


class TestRepartitionForJoin:
    """Test the repartition_for_join function."""

    def test_returns_list_of_intervals(self, int_keyed_ht):
        """Test that the function returns a non-empty list of Intervals."""
        intervals = repartition_for_join(int_keyed_ht, n_partitions=5)

        assert isinstance(intervals, list)
        assert len(intervals) > 0
        assert all(isinstance(i, Interval) for i in intervals)

    def test_n_partitions_controls_count(self, int_keyed_ht):
        """Test that `n_partitions` sets the number of calculated intervals."""
        intervals = repartition_for_join(int_keyed_ht, n_partitions=5)

        assert len(intervals) == 5

    def test_n_partitions_overrides_percent(self, int_keyed_ht):
        """Test that `n_partitions` takes precedence over `new_partition_percent`."""
        intervals = repartition_for_join(
            int_keyed_ht, new_partition_percent=5.0, n_partitions=3
        )

        # The count should follow n_partitions (3), not the percent-derived value
        # (int(10 * 5.0) == 50).
        assert len(intervals) == 3

    def test_new_partition_percent_default(self, int_keyed_ht):
        """Test the interval count derived from the default `new_partition_percent`."""
        intervals = repartition_for_join(int_keyed_ht)

        # The default requests int(ht.n_partitions() * 1.1) == int(10 * 1.1) == 11
        # partitions; the partitioner may return slightly fewer (here, 10).
        assert len(intervals) == 10

    def test_new_partition_percent_custom(self, int_keyed_ht):
        """Test the interval count derived from a custom `new_partition_percent`."""
        intervals = repartition_for_join(int_keyed_ht, new_partition_percent=2.0)

        # int(ht.n_partitions() * 2.0) == int(10 * 2.0) == 20.
        assert len(intervals) == 20

    def test_percent_less_than_one_warns(self, int_keyed_ht, caplog):
        """Test that a `new_partition_percent` below 1 logs a warning."""
        with caplog.at_level(logging.WARNING, logger="gnomad.utils.file_utils"):
            repartition_for_join(int_keyed_ht, new_partition_percent=0.5)

        assert any(
            "new_partition_percent value is less than 1" in record.message
            for record in caplog.records
        )

    def test_reads_table_from_path(self, int_keyed_ht, tmp_path):
        """Test that a string path is read with `hl.read_table`."""
        path = str(tmp_path / "int_keyed.ht")
        int_keyed_ht.write(path, overwrite=True)

        intervals = repartition_for_join(path, n_partitions=5)

        assert len(intervals) == 5
        assert all(isinstance(i, Interval) for i in intervals)

    def test_default_intervals_are_key_structs(self, locus_keyed_ht):
        """Test that key-struct intervals are returned by default."""
        intervals = repartition_for_join(locus_keyed_ht, n_partitions=5)

        assert all(isinstance(i.start, Struct) for i in intervals)
        assert all("locus" in i.start for i in intervals)

    def test_locus_intervals(self, locus_keyed_ht):
        """Test that `locus_intervals=True` returns bare-Locus intervals."""
        intervals = repartition_for_join(
            locus_keyed_ht, n_partitions=5, locus_intervals=True
        )

        assert all(isinstance(i.start, Locus) for i in intervals)
        assert all(isinstance(i.end, Locus) for i in intervals)

    def test_locus_intervals_requires_locus_key(self, int_keyed_ht):
        """Test that `locus_intervals=True` needs a locus-keyed Table."""
        # The docstring requires the input to be keyed by `locus`; an `idx`-keyed
        # Table has no `locus` field in its interval bounds.
        with pytest.raises(KeyError):
            repartition_for_join(int_keyed_ht, n_partitions=5, locus_intervals=True)

    def test_intervals_usable_for_copartitioned_read(self, int_keyed_ht, tmp_path):
        """Test that the returned intervals work as co-partition read intervals.

        This is the function's actual purpose: the intervals are passed as
        `_intervals` when reading a Table so that joins are co-partitioned.
        """
        path = str(tmp_path / "roundtrip.ht")
        int_keyed_ht.write(path, overwrite=True)

        intervals = repartition_for_join(path, n_partitions=5)
        repartitioned = hl.read_table(path, _intervals=intervals)

        # Reading with the intervals preserves all rows and yields exactly the
        # requested number of partitions.
        assert repartitioned.count() == int_keyed_ht.count()
        assert repartitioned.n_partitions() == 5
