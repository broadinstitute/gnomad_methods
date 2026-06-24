"""Tests for the file_utils utility module."""

import os
from unittest.mock import patch

import hail as hl

from gnomad.utils.file_utils import file_exists, read_list_data


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

    @patch("gnomad.utils.file_utils.hfs.exists")
    def test_gcs_hail_table_checks_success_path(self, mock_exists):
        """Test the `_SUCCESS` path built for a gs:// Hail Table (no network)."""
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
