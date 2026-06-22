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


@pytest.fixture(scope="module")
def locus_alleles_keyed_ht():
    """Return a Table keyed by `locus` and `alleles` with 10 partitions."""
    ht = hl.utils.range_table(100, n_partitions=10)
    ht = ht.annotate(
        locus=hl.locus("chr1", ht.idx + 1, reference_genome="GRCh38"),
        alleles=hl.array(["A", "C"]),
    )
    return ht.key_by("locus", "alleles")


@pytest.fixture(scope="module")
def synthetic_vds():
    """Build a tiny VDS (20 sites x 4 samples) for interval-based read tests."""
    samples = [f"s{i}" for i in range(4)]
    positions = list(range(1000, 1000 + 20 * 100, 100))

    vd = hl.Table.parallelize(
        [
            {
                "locus": hl.locus("chr1", p, reference_genome="GRCh38"),
                "alleles": ["A", "T"],
                "s": s,
                "GT": hl.call(0, 1),
            }
            for p in positions
            for s in samples
        ],
        hl.tstruct(
            locus=hl.tlocus("GRCh38"),
            alleles=hl.tarray(hl.tstr),
            s=hl.tstr,
            GT=hl.tcall,
        ),
    ).to_matrix_table(row_key=["locus", "alleles"], col_key=["s"])

    # reference_data: 1-base ref blocks at every site (END == position).
    rd = hl.Table.parallelize(
        [
            {
                "locus": hl.locus("chr1", p, reference_genome="GRCh38"),
                "s": s,
                "END": p,
                "LEN": 1,
            }
            for p in positions
            for s in samples
        ],
        hl.tstruct(locus=hl.tlocus("GRCh38"), s=hl.tstr, END=hl.tint32, LEN=hl.tint32),
    ).to_matrix_table(row_key=["locus"], col_key=["s"])

    return hl.vds.VariantDataset(reference_data=rd, variant_data=vd)


class TestRepartitionForJoin:
    """Test the repartition_for_join function."""

    def test_returns_list_of_intervals(self, int_keyed_ht):
        """Test that the function returns a non-empty list of Intervals."""
        intervals = repartition_for_join(int_keyed_ht, n_partitions=5)

        assert isinstance(intervals, list)
        assert all(isinstance(i, Interval) for i in intervals)

    def test_intervals_are_sorted_and_non_overlapping(self, int_keyed_ht):
        """Test that the intervals are ordered, contiguous, and non-overlapping.

        Full coverage of the key space is verified separately by the round-trip
        read tests (which confirm no rows are dropped).
        """
        intervals = repartition_for_join(int_keyed_ht, n_partitions=5)

        starts = [i.start["idx"] for i in intervals]
        ends = [i.end["idx"] for i in intervals]

        # Starts are strictly increasing and each interval ends where the next
        # begins, so the intervals are sorted, non-overlapping, and contiguous.
        assert starts == sorted(starts)
        assert all(ends[k] == starts[k + 1] for k in range(len(intervals) - 1))

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

    # None passes through the default (1.1); a value exercises a custom percent.
    @pytest.mark.parametrize("percent", [None, 2.0])
    def test_new_partition_percent_controls_count(self, int_keyed_ht, percent):
        """Test the interval count derived from `new_partition_percent`."""
        orig = int_keyed_ht.n_partitions()
        kwargs = {} if percent is None else {"new_partition_percent": percent}
        intervals = repartition_for_join(int_keyed_ht, **kwargs)

        # `_calculate_new_partitions` is not guaranteed to return exactly the
        # requested count (int(orig * percent)). With a percent > 1 the result is
        # between the original partition count and the requested count.
        requested = int(orig * (percent if percent is not None else 1.1))
        assert orig <= len(intervals) <= requested

    def test_percent_less_than_one_warns(self, int_keyed_ht, caplog):
        """Test that a `new_partition_percent` below 1 logs a warning."""
        with caplog.at_level(logging.WARNING, logger="gnomad.utils.file_utils"):
            repartition_for_join(int_keyed_ht, new_partition_percent=0.5)

        assert any(
            "new_partition_percent value is less than 1" in record.message
            for record in caplog.records
        )

    def test_path_and_table_inputs_match(self, int_keyed_ht, tmp_path):
        """Test that a path string and an in-memory Table give the same intervals."""
        path = str(tmp_path / "int_keyed.ht")
        int_keyed_ht.write(path, overwrite=True)

        from_path = repartition_for_join(path, n_partitions=5)
        from_table = repartition_for_join(int_keyed_ht, n_partitions=5)

        assert all(isinstance(i, Interval) for i in from_path)
        assert from_path == from_table

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

    def test_locus_alleles_default_intervals(self, locus_alleles_keyed_ht):
        """Test key-struct intervals over a `locus`/`alleles` keyed Table."""
        intervals = repartition_for_join(locus_alleles_keyed_ht, n_partitions=5)

        assert all(isinstance(i.start, Struct) for i in intervals)
        assert all("locus" in i.start and "alleles" in i.start for i in intervals)

    def test_locus_alleles_locus_intervals(self, locus_alleles_keyed_ht):
        """Test `locus_intervals=True` over a `locus`/`alleles` keyed Table."""
        intervals = repartition_for_join(
            locus_alleles_keyed_ht, n_partitions=5, locus_intervals=True
        )

        # Only the `locus` component is kept, as bare Locus interval bounds.
        assert all(isinstance(i.start, Locus) for i in intervals)
        assert all(isinstance(i.end, Locus) for i in intervals)

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

    def test_intervals_usable_for_matrixtable_read(self, tmp_path):
        """Test that the intervals work as co-partition read intervals for an MT.

        Intervals are computed from the MatrixTable's `rows()` Table and passed
        to `hl.read_matrix_table(..., _intervals=...)`.
        """
        path = str(tmp_path / "roundtrip.mt")
        mt = hl.utils.range_matrix_table(100, 10, n_partitions=10)
        mt = mt.key_rows_by(
            locus=hl.locus("chr1", mt.row_idx + 1, reference_genome="GRCh38")
        )
        mt.write(path, overwrite=True)
        mt = hl.read_matrix_table(path)

        intervals = repartition_for_join(mt.rows(), n_partitions=5)
        repartitioned = hl.read_matrix_table(path, _intervals=intervals)

        assert repartitioned.count_rows() == mt.count_rows()
        assert repartitioned.n_partitions() == 5

    def test_intervals_usable_for_vds_read(self, synthetic_vds, tmp_path):
        """Test that locus intervals work as co-partition read intervals for a VDS.

        Locus intervals are computed from the VDS `variant_data.rows()` Table and
        passed to `hl.vds.read_vds(..., intervals=...)`, the VDS reader's expected
        form (hence `locus_intervals=True`).
        """
        path = str(tmp_path / "roundtrip.vds")
        synthetic_vds.write(path, overwrite=True)
        vds = hl.vds.read_vds(path)

        intervals = repartition_for_join(
            vds.variant_data.rows(), n_partitions=5, locus_intervals=True
        )
        repartitioned = hl.vds.read_vds(path, intervals=intervals)

        assert (
            repartitioned.variant_data.count_rows()
            == synthetic_vds.variant_data.count_rows()
        )
        assert repartitioned.variant_data.n_partitions() == 5
