"""Tests for resource classes."""

from typing import List, Tuple, Union
from unittest.mock import patch

import pytest

from gnomad.resources import resource_utils
from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)


class TestTableResource:
    """Tests for TableResource."""

    @patch("hail.read_table")
    def test_read_table(self, read_table):
        """Test that Table is read from path."""
        resource = resource_utils.TableResource("gs://gnomad-public/table.ht")

        ds = resource.ht()
        read_table.assert_called_with("gs://gnomad-public/table.ht")
        assert ds == read_table.return_value


class TestMatrixTableResource:
    """Tests for MatrixTableResource."""

    @patch("hail.read_matrix_table")
    def test_read_matrix_table(self, read_matrix_table):
        """Test that MatrixTable is read from path."""
        resource = resource_utils.MatrixTableResource(
            "gs://gnomad-public/matrix_table.mt"
        )

        ds = resource.mt()
        read_matrix_table.assert_called_with("gs://gnomad-public/matrix_table.mt")
        assert ds == read_matrix_table.return_value


class TestPedigreeResource:
    """Tests for PedigreeResource."""

    @patch("hail.Pedigree.read")
    def test_read_pedigree(self, read_pedigree):
        """Test that Pedigree is read from path."""
        resource = resource_utils.PedigreeResource("gs://gnomad-public/pedigree.ped")

        ds = resource.pedigree()
        read_pedigree.assert_called()
        print(read_pedigree.call_args)
        assert read_pedigree.call_args[0][0] == "gs://gnomad-public/pedigree.ped"
        assert ds == read_pedigree.return_value

    @patch("hail.import_fam")
    def test_read_fam(self, import_fam):
        """Test that Table is imported from path."""
        resource = resource_utils.PedigreeResource("gs://gnomad-public/pedigree.fam")

        ds = resource.ht()
        import_fam.assert_called()
        assert import_fam.call_args[0][0] == "gs://gnomad-public/pedigree.fam"
        assert ds == import_fam.return_value


class TestBlockMatrixResource:
    """Tests for BlockMatrixResource."""

    @patch("hail.linalg.BlockMatrix.read")
    def test_read_block_matrix(self, read_block_matrix):
        """Test that BlockMatrix is read from path."""
        resource = resource_utils.BlockMatrixResource(
            "gs://gnomad-public/block_matrix.bm"
        )

        ds = resource.bm()
        read_block_matrix.assert_called_with("gs://gnomad-public/block_matrix.bm")
        assert ds == read_block_matrix.return_value


def gnomad_public_resource_test_parameters(
    path: str,
) -> List[Tuple[str, Union[GnomadPublicResourceSource, str], str]]:
    """
    Get parameters for gnomAD public resource tests.

    :param path: Path to resource file inside gnomAD bucket.
    """
    return [
        (
            f"gs://gnomad-public{path}",
            GnomadPublicResourceSource.GNOMAD,
            f"gs://gnomad-public{path}",
        ),
        (
            f"gs://gnomad-public-requester-pays{path}",
            GnomadPublicResourceSource.GNOMAD,
            f"gs://gnomad-public-requester-pays{path}",
        ),
        (
            f"gs://gnomad-public{path}",
            GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS,
            f"gs://gcp-public-data--gnomad{path}",
        ),
        (
            f"gs://gnomad-public-requester-pays{path}",
            GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS,
            f"gs://gcp-public-data--gnomad{path}",
        ),
        (
            f"gs://gnomad-public{path}",
            "gs://my-bucket/gnomad-resources",
            f"gs://my-bucket/gnomad-resources{path}",
        ),
        (
            f"gs://gnomad-public-requester-pays{path}",
            "gs://my-bucket/gnomad-resources",
            f"gs://my-bucket/gnomad-resources{path}",
        ),
    ]


class TestGnomadPublicTableResource:
    """Tests for GnomadPublicTableResource."""

    @pytest.mark.parametrize(
        "resource_path,source,expected_read_path",
        gnomad_public_resource_test_parameters("/table.ht"),
    )
    @patch("hail.read_table")
    def test_read_gnomad_public_table_resource(
        self, read_table, resource_path, source, expected_read_path
    ):
        """Test that Table can be read from different sources."""
        resource = resource_utils.GnomadPublicTableResource(resource_path)

        gnomad_public_resource_configuration.source = source

        resource.ht()
        read_table.assert_called_with(expected_read_path)


class TestGnomadPublicMatrixTableResource:
    """Tests for GnomadPublicMatrixTableResource."""

    @pytest.mark.parametrize(
        "resource_path,source,expected_read_path",
        gnomad_public_resource_test_parameters("/matrix_table.mt"),
    )
    @patch("hail.read_matrix_table")
    def test_read_gnomad_public_matrix_table_resource(
        self, read_matrix_table, resource_path, source, expected_read_path
    ):
        """Test that MatrixTable can be read from different sources."""
        resource = resource_utils.GnomadPublicMatrixTableResource(resource_path)

        gnomad_public_resource_configuration.source = source

        resource.mt()
        read_matrix_table.assert_called_with(expected_read_path)


class TestGnomadPublicPedigreeResource:
    """Tests for GnomadPublicPedigreeResource."""

    @pytest.mark.parametrize(
        "resource_path,source,expected_read_path",
        gnomad_public_resource_test_parameters("/pedigree.ped"),
    )
    @patch("hail.Pedigree.read")
    def test_read_gnomad_public_pedigree_resource(
        self, read_pedigree, resource_path, source, expected_read_path
    ):
        """Test that Pedigree can be read from different sources."""
        resource = resource_utils.GnomadPublicPedigreeResource(resource_path)

        gnomad_public_resource_configuration.source = source

        resource.pedigree()
        read_pedigree.assert_called()
        assert read_pedigree.call_args[0][0] == expected_read_path

    @pytest.mark.parametrize(
        "resource_path,source,expected_read_path",
        gnomad_public_resource_test_parameters("/pedigree.fam"),
    )
    @patch("hail.import_fam")
    def test_import_gnomad_public_pedigree_resource(
        self, import_fam, resource_path, source, expected_read_path
    ):
        """Test that pedigree can be imported from different sources."""
        resource = resource_utils.GnomadPublicPedigreeResource(resource_path)

        gnomad_public_resource_configuration.source = source

        resource.ht()
        import_fam.assert_called()
        assert import_fam.call_args[0][0] == expected_read_path


class TestGnomadPublicBlockMatrixResource:
    """Tests for GnomadPublicBlockMatrixResource."""

    @pytest.mark.parametrize(
        "resource_path,source,expected_read_path",
        gnomad_public_resource_test_parameters("/block_matrix.bm"),
    )
    @patch("hail.linalg.BlockMatrix.read")
    def test_read_gnomad_public_block_matrix_resource(
        self, read_block_matrix, resource_path, source, expected_read_path
    ):
        """Test that BlockMatrix can be read from different sources."""
        resource = resource_utils.GnomadPublicBlockMatrixResource(resource_path)

        gnomad_public_resource_configuration.source = source

        resource.bm()
        read_block_matrix.assert_called_with(expected_read_path)
