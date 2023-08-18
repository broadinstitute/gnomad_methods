"""Tests for resource classes."""

import os
from typing import List, Tuple, Union
from unittest.mock import patch

import pytest

from gnomad.resources import resource_utils
from gnomad.resources.config import (
    GnomadPublicResourceSource,
    get_default_public_resource_source,
    gnomad_public_resource_configuration,
)


class TestTableResource:
    """Tests for TableResource."""

    @patch("hail.read_table")
    def test_read_table(self, read_table):
        """Test that Table is read from path."""
        resource = resource_utils.TableResource(
            "gs://gnomad-public-requester-pays/table.ht"
        )

        ds = resource.ht()
        read_table.assert_called_with("gs://gnomad-public-requester-pays/table.ht")
        assert ds == read_table.return_value


class TestMatrixTableResource:
    """Tests for MatrixTableResource."""

    @patch("hail.read_matrix_table")
    def test_read_matrix_table(self, read_matrix_table):
        """Test that MatrixTable is read from path."""
        resource = resource_utils.MatrixTableResource(
            "gs://gnomad-public-requester-pays/matrix_table.mt"
        )

        ds = resource.mt()
        read_matrix_table.assert_called_with(
            "gs://gnomad-public-requester-pays/matrix_table.mt"
        )
        assert ds == read_matrix_table.return_value


class TestPedigreeResource:
    """Tests for PedigreeResource."""

    @patch("hail.Pedigree.read")
    def test_read_pedigree(self, read_pedigree):
        """Test that Pedigree is read from path."""
        resource = resource_utils.PedigreeResource(
            "gs://gnomad-public-requester-pays/pedigree.ped"
        )

        ds = resource.pedigree()
        read_pedigree.assert_called()
        print(read_pedigree.call_args)
        assert (
            read_pedigree.call_args[0][0]
            == "gs://gnomad-public-requester-pays/pedigree.ped"
        )
        assert ds == read_pedigree.return_value

    @patch("hail.import_fam")
    def test_read_fam(self, import_fam):
        """Test that Table is imported from path."""
        resource = resource_utils.PedigreeResource(
            "gs://gnomad-public-requester-pays/pedigree.fam"
        )

        ds = resource.ht()
        import_fam.assert_called()
        assert (
            import_fam.call_args[0][0]
            == "gs://gnomad-public-requester-pays/pedigree.fam"
        )
        assert ds == import_fam.return_value


class TestBlockMatrixResource:
    """Tests for BlockMatrixResource."""

    @patch("hail.linalg.BlockMatrix.read")
    def test_read_block_matrix(self, read_block_matrix):
        """Test that BlockMatrix is read from path."""
        resource = resource_utils.BlockMatrixResource(
            "gs://gnomad-public-requester-pays/block_matrix.bm"
        )

        ds = resource.bm()
        read_block_matrix.assert_called_with(
            "gs://gnomad-public-requester-pays/block_matrix.bm"
        )
        assert ds == read_block_matrix.return_value


class TestDefaultPublicResourceSource:
    """Tests for default public resource source."""

    @pytest.mark.parametrize(
        "default_source,expected_path",
        [
            (
                GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS,
                "gs://gcp-public-data--gnomad/example.ht",
            ),
            (
                GnomadPublicResourceSource.REGISTRY_OF_OPEN_DATA_ON_AWS,
                "s3a://gnomad-public-us-east-1/example.ht",
            ),
            (
                GnomadPublicResourceSource.AZURE_OPEN_DATASETS,
                "wasbs://dataset@datasetgnomad.blob.core.windows.net/example.ht",
            ),
            (
                "gs://my-bucket/gnomad-resources",
                "gs://my-bucket/gnomad-resources/example.ht",
            ),
        ],
    )
    def test_read_from_default_source(self, default_source, expected_path):
        """Test that resource paths use default source when no source is configured."""
        gnomad_public_resource_configuration._source = None

        with patch(
            "gnomad.resources.config.get_default_public_resource_source",
            return_value=default_source,
        ):
            resource = resource_utils.GnomadPublicTableResource(
                "gs://gnomad-public-requester-pays/example.ht"
            )
            assert resource.path == expected_path

    @pytest.mark.parametrize(
        "configured_default_source,expected_default_source",
        [
            ("gnomAD", GnomadPublicResourceSource.GNOMAD),
            (
                "Google Cloud Public Datasets",
                GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS,
            ),
            (
                "Registry of Open Data on AWS",
                GnomadPublicResourceSource.REGISTRY_OF_OPEN_DATA_ON_AWS,
            ),
            ("Azure Open Datasets", GnomadPublicResourceSource.AZURE_OPEN_DATASETS),
            ("gs://my-bucket/gnomad-resources", "gs://my-bucket/gnomad-resources"),
        ],
    )
    def test_get_default_source_from_environment(
        self, configured_default_source, expected_default_source
    ):
        """Test that default source is read from environment variable."""
        with patch.dict(
            os.environ,
            {"GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE": configured_default_source},
        ):
            assert get_default_public_resource_source() == expected_default_source

    @pytest.mark.parametrize(
        "cloud_spark_provider,expected_default_source",
        [
            ("dataproc", GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS),
            ("hdinsight", GnomadPublicResourceSource.AZURE_OPEN_DATASETS),
            ("unknown", GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS),
            (None, GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS),
        ],
    )
    def test_get_default_source_from_cloud_spark_provider(
        self, cloud_spark_provider, expected_default_source
    ):
        """Test that default source is set based on cloud Spark provider."""
        with patch(
            "hail.utils.guess_cloud_spark_provider",
            return_value=cloud_spark_provider,
            create=True,
        ):
            assert get_default_public_resource_source() == expected_default_source

    def test_default_source_from_environment_overrides_cloud_spark_provider(self):
        """
        Test that a default source configured in environment variables is preferred.

        Make sure the environment variables is used over the one for the current cloud
        Spark provider.
        """
        with (
            patch(
                "hail.utils.guess_cloud_spark_provider",
                return_value="hdinsight",
                create=True,
            ),
            patch.dict(
                os.environ,
                {
                    "GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE": (
                        "gs://my-bucket/gnomad-resources"
                    )
                },
            ),
        ):
            assert (
                get_default_public_resource_source()
                == "gs://my-bucket/gnomad-resources"
            )


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
            GnomadPublicResourceSource.REGISTRY_OF_OPEN_DATA_ON_AWS,
            f"s3a://gnomad-public-us-east-1{path}",
        ),
        (
            f"gs://gnomad-public-requester-pays{path}",
            GnomadPublicResourceSource.REGISTRY_OF_OPEN_DATA_ON_AWS,
            f"s3a://gnomad-public-us-east-1{path}",
        ),
        (
            f"gs://gnomad-public{path}",
            GnomadPublicResourceSource.AZURE_OPEN_DATASETS,
            f"wasbs://dataset@datasetgnomad.blob.core.windows.net{path}",
        ),
        (
            f"gs://gnomad-public-requester-pays{path}",
            GnomadPublicResourceSource.AZURE_OPEN_DATASETS,
            f"wasbs://dataset@datasetgnomad.blob.core.windows.net{path}",
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

        with patch.object(resource, "is_resource_available", return_value=True):
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

        with patch.object(resource, "is_resource_available", return_value=True):
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

        with patch.object(resource, "is_resource_available", return_value=True):
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

        with patch.object(resource, "is_resource_available", return_value=True):
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

        with patch.object(resource, "is_resource_available", return_value=True):
            resource.bm()
            read_block_matrix.assert_called_with(expected_read_path)
