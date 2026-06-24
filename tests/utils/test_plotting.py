"""Tests for the plotting utility module."""

import os

import hail as hl

from gnomad.utils.plotting import _ls, get_rows_data


class TestLs:
    """Test the _ls helper, which mirrors the deprecated hl.hadoop_ls."""

    def test_excludes_hidden_dotfiles(self, tmp_path):
        """Test that dot-prefixed files (e.g. local `.crc`) are excluded.

        `hfs.ls` lists hidden files that `hl.hadoop_ls` did not; writing a Hail
        Table on the local filesystem creates `.crc` checksum files, which must
        not leak through `_ls`.
        """
        path = str(tmp_path / "t.ht")
        hl.utils.range_table(50, n_partitions=3).write(path)
        parts = os.path.join(path, "rows", "parts")

        names = [os.path.basename(e["path"]) for e in _ls(parts)]

        assert names, "expected part files to be listed"
        assert all(not n.startswith(".") for n in names)
        assert not any(n.endswith(".crc") for n in names)

    def test_keeps_underscore_success(self, tmp_path):
        """Test that `_SUCCESS` (underscore-prefixed) is retained, like hl.hadoop_ls."""
        path = str(tmp_path / "t.ht")
        hl.utils.range_table(50, n_partitions=3).write(path)

        names = [os.path.basename(e["path"]) for e in _ls(path)]

        assert "_SUCCESS" in names


class TestGetRowsData:
    """Test get_rows_data partition-file parsing."""

    def test_returns_one_size_per_partition(self, tmp_path):
        """Test that part files are parsed without tripping on hidden `.crc` files.

        Regression guard: when `ls` includes `.crc` files, the `part-` parsing
        raised IndexError and file-size counts were inflated.
        """
        path = str(tmp_path / "t.ht")
        hl.utils.range_table(100, n_partitions=5).write(path)
        rows_files = _ls(os.path.join(path, "rows"))

        _, file_sizes = get_rows_data(rows_files)

        # Exactly one size per partition; the `.crc` files were not counted.
        assert len(file_sizes) == 5
        assert all(s > 0 for s in file_sizes)
