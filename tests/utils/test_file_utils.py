"""Tests for the file_utils utility module."""

import gzip
import os
from unittest.mock import patch

import hail as hl
import pytest

from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import (
    check_file_exists_raise_error,
    file_exists,
    read_list_data,
)


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
