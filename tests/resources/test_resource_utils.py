"""Tests for resource classes."""

from unittest.mock import patch

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


class TestGnomadPublicTableResource:
    """Tests for GnomadPublicTableResource."""

    @patch("hail.read_table")
    def test_gnomad_public_table_resource(self, read_table):
        """Test that Table can be read from gnomAD bucket."""
        resource = resource_utils.GnomadPublicTableResource(
            "gs://gnomad-public/table.ht"
        )

        gnomad_public_resource_configuration.source = GnomadPublicResourceSource.GNOMAD

        resource.ht()
        read_table.assert_called_with("gs://gnomad-public/table.ht")

    @patch("hail.read_table")
    def test_gnomad_public_table_resource_google_cloud_public_datasets(
        self, read_table
    ):
        """Test that Table can be read from Google Cloud Public Datasets bucket."""
        resource = resource_utils.GnomadPublicTableResource(
            "gs://gnomad-public/table.ht"
        )

        gnomad_public_resource_configuration.source = (
            GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
        )

        resource.ht()
        read_table.assert_called_with("gs://gcp-public-data--gnomad/table.ht")

    @patch("hail.read_table")
    def test_gnomad_public_table_resource_custom_source(self, read_table):
        """Test that Table can be read from custom source."""
        resource = resource_utils.GnomadPublicTableResource(
            "gs://gnomad-public/table.ht"
        )

        gnomad_public_resource_configuration.source = "gs://my-bucket/gnomad-resources"

        resource.ht()
        read_table.assert_called_with("gs://my-bucket/gnomad-resources/table.ht")
