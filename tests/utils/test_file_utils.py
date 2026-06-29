"""Tests for the file_utils utility module."""

import gzip
import logging
import os
from unittest.mock import patch

import hail as hl
import pytest
from hail.genetics import Locus
from hail.utils import Interval, Struct

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import (
    check_file_exists_raise_error,
    file_exists,
    read_list_data,
    repartition_for_join,
)

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


class TestFileExists:
    """Test the file_exists function."""

    def test_plain_file_exists(self, tmp_path):
        """Test that an existing plain file is detected."""
        path = str(tmp_path / "data.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("content")

        assert file_exists(path) is True

    def test_missing_file(self, tmp_path):
        """Test that a missing path is reported as absent."""
        assert file_exists(str(tmp_path / "missing.txt")) is False

    def test_hail_table_checks_success_file(self, tmp_path):
        """Test that a written `.ht` is detected via its `_SUCCESS` file."""
        path = str(tmp_path / "t.ht")
        hl.utils.range_table(10).write(path)

        assert file_exists(path) is True

    def test_hail_table_missing_success_file(self, tmp_path):
        """Test that an `.ht` directory without `_SUCCESS` is not detected."""
        path = str(tmp_path / "t.ht")
        hl.utils.range_table(10).write(path)
        os.remove(os.path.join(path, "_SUCCESS"))

        assert file_exists(path) is False

    # No network needed for mock gs:// paths
    @patch("gnomad.utils.file_utils.hfs.exists")
    def test_gcs_hail_table_checks_success_path(self, mock_exists):
        """Test the `_SUCCESS` path built for a gs:// Hail Table."""
        mock_exists.return_value = True

        assert file_exists("gs://my-bucket/t.ht") is True
        mock_exists.assert_called_once_with("gs://my-bucket/t.ht/_SUCCESS")

    @patch("gnomad.utils.file_utils.hfs.exists")
    def test_gcs_vds_checks_both_success_paths(self, mock_exists):
        """Test that a gs:// VDS checks both component `_SUCCESS` files."""
        mock_exists.return_value = True

        assert file_exists("gs://my-bucket/data.vds") is True
        assert mock_exists.call_count == 2
        mock_exists.assert_any_call("gs://my-bucket/data.vds/reference_data/_SUCCESS")
        mock_exists.assert_any_call("gs://my-bucket/data.vds/variant_data/_SUCCESS")

    @patch("gnomad.utils.file_utils.hfs.exists")
    def test_gcs_vds_missing_one_component_is_false(self, mock_exists):
        """Test that a VDS missing one component `_SUCCESS` is reported absent."""
        # reference_data present, variant_data missing -> overall False via all().
        mock_exists.side_effect = lambda p: p.endswith("reference_data/_SUCCESS")

        assert file_exists("gs://my-bucket/data.vds") is False

    @patch("gnomad.utils.file_utils.hfs.exists")
    def test_gcs_plain_path_checked_directly(self, mock_exists):
        """Test that a gs:// non-Hail path is checked as-is (no `_SUCCESS`)."""
        mock_exists.return_value = False

        assert file_exists("gs://my-bucket/list.txt") is False
        mock_exists.assert_called_once_with("gs://my-bucket/list.txt")


class TestReadListData:
    """Test the read_list_data function."""

    def test_reads_lines_stripped(self, tmp_path):
        """Test that each line is read into a list with whitespace stripped."""
        path = str(tmp_path / "list.txt")
        with open(path, "w", encoding="utf-8") as f:
            f.write("a\nb \n c\n")

        assert read_list_data(path) == ["a", "b", "c"]

    def test_reads_gzipped_lines_stripped(self, tmp_path):
        """Test that a .gz file is read in text mode with whitespace stripped."""
        path = str(tmp_path / "list.txt.gz")
        with gzip.open(path, mode="wt", encoding="utf-8") as f:
            f.write("a\nb \n c\n")

        assert read_list_data(path) == ["a", "b", "c"]


class TestCheckFileExistsRaiseError:
    """Test the check_file_exists_raise_error function."""

    @patch("gnomad.utils.file_utils.file_exists")
    def test_returns_all_exist_without_error_flags(self, mock_exists):
        """Test the all-exist boolean is returned (and a str path is one file)."""
        mock_exists.return_value = True
        assert check_file_exists_raise_error("a.ht") is True

        mock_exists.return_value = False
        assert check_file_exists_raise_error(["a.ht", "b.ht"]) is False

    @patch("gnomad.utils.file_utils.file_exists")
    def test_raises_when_existing_file_and_error_if_exists(self, mock_exists):
        """Test that an existing file raises, naming only the offending file."""
        mock_exists.side_effect = lambda f: f == "exists.ht"
        with pytest.raises(DataException, match="already exist") as exc_info:
            check_file_exists_raise_error(
                ["exists.ht", "missing.ht"], error_if_exists=True
            )
        assert "exists.ht" in str(exc_info.value)
        assert "missing.ht" not in str(exc_info.value)

    @patch("gnomad.utils.file_utils.file_exists")
    def test_raises_when_missing_file_and_error_if_not_exists(self, mock_exists):
        """Test that a missing file raises when error_if_not_exists is set."""
        mock_exists.return_value = False
        with pytest.raises(DataException, match="do not exist.*a.ht"):
            check_file_exists_raise_error("a.ht", error_if_not_exists=True)

    @patch("gnomad.utils.file_utils.file_exists")
    def test_raises_with_both_messages_when_both_flags_set(self, mock_exists):
        """Test that both flags together report existing and missing files in one error."""
        mock_exists.side_effect = lambda f: f == "exists.ht"
        with pytest.raises(DataException) as exc_info:
            check_file_exists_raise_error(
                ["exists.ht", "missing.ht"],
                error_if_exists=True,
                error_if_not_exists=True,
            )
        msg = str(exc_info.value)
        # Both halves of the combined message are present, each naming its file.
        assert "already exist" in msg and "exists.ht" in msg
        assert "do not exist" in msg and "missing.ht" in msg

    @patch("gnomad.utils.file_utils.file_exists")
    def test_no_raise_when_flag_set_but_condition_not_met(self, mock_exists):
        """Test that a set flag only raises when its condition is actually triggered."""
        # error_if_exists is set but nothing exists: no raise, returns False.
        mock_exists.return_value = False
        assert (
            check_file_exists_raise_error(["a.ht", "b.ht"], error_if_exists=True)
            is False
        )
        # error_if_not_exists is set but everything exists: no raise, returns True.
        mock_exists.return_value = True
        assert (
            check_file_exists_raise_error(["a.ht", "b.ht"], error_if_not_exists=True)
            is True
        )

    @patch("gnomad.utils.file_utils.file_exists")
    def test_custom_error_messages_are_used(self, mock_exists):
        """Test that caller-supplied error message prefixes are honored."""
        mock_exists.return_value = True
        with pytest.raises(DataException, match="CUSTOM EXISTS: a.ht"):
            check_file_exists_raise_error(
                "a.ht",
                error_if_exists=True,
                error_if_exists_msg="CUSTOM EXISTS: ",
            )
