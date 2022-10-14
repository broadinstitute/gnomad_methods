"""Configuration for loading resources."""

import logging
import os
from enum import Enum
from typing import Union

logger = logging.getLogger(__name__)


class GnomadPublicResourceSource(Enum):
    """Sources for public gnomAD resources."""

    GNOMAD = "gnomAD"
    GOOGLE_CLOUD_PUBLIC_DATASETS = "Google Cloud Public Datasets"
    REGISTRY_OF_OPEN_DATA_ON_AWS = "Registry of Open Data on AWS"
    AZURE_OPEN_DATASETS = "Azure Open Datasets"


def get_default_public_resource_source() -> Union[GnomadPublicResourceSource, str]:
    """
    Get the default source for public gnomAD resources.

    The default source is determined by...

    - If the ``GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE`` environment variable is set, use the source configured there.
    - Otherwise, if Hail determines that is is running in a cloud provider's Spark environment, use the source from that cloud provider.
      For example, use Azure Open Datasets if running on an Azure HDInsight cluster.
    - Otherwise, use Google Cloud Public Datasets.

    :returns: Default resource source
    """
    default_source_from_env = os.getenv("GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE", None)
    if default_source_from_env:
        # Convert to a GnomadPublicResourceSource enum if possible
        try:
            default_source = GnomadPublicResourceSource(default_source_from_env)
            logger.info(
                "Using configured source for gnomAD resources: %s", default_source.value
            )
            return default_source
        except ValueError:
            logger.info(
                "Using configured custom source for gnomAD resources: %s",
                default_source_from_env,
            )
            return default_source_from_env

    try:
        from hail.utils import guess_cloud_spark_provider
    except ImportError:
        pass
    else:
        cloud_spark_provider = guess_cloud_spark_provider()
        default_resource_sources_by_provider = {
            "dataproc": GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS,
            "hdinsight": GnomadPublicResourceSource.AZURE_OPEN_DATASETS,
        }
        if cloud_spark_provider:
            try:
                default_source_from_provider = default_resource_sources_by_provider[
                    cloud_spark_provider
                ]
                logger.info(
                    "Using default source for gnomAD resources based on cloud"
                    " provider: %s",
                    default_source_from_provider,
                )
                return default_source_from_provider
            except KeyError:
                pass

    return GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS


class _GnomadPublicResourceConfiguration:
    """Configuration for public gnomAD resources."""

    _source: Union[GnomadPublicResourceSource, str, None] = None

    @property
    def source(self) -> Union[GnomadPublicResourceSource, str]:
        """
        Get the source for public gnomAD resource files.

        This is used to determine which URLs gnomAD resources will be loaded from.

        :returns: Source name or path to root of resources directory
        """
        if self._source is None:
            self._source = get_default_public_resource_source()

        return self._source

    @source.setter
    def source(self, source: Union[GnomadPublicResourceSource, str]) -> None:
        """
        Set the default source for resource files.

        This is used to determine which URLs gnomAD resources will be loaded from.

        :param source: Source name or path to root of resources directory
        """
        self._source = source


gnomad_public_resource_configuration = _GnomadPublicResourceConfiguration()
