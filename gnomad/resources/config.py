"""Configuration for loading resources."""

import os
from enum import Enum
from typing import Union


class GnomadPublicResourceSource(Enum):
    """Sources for public gnomAD resources."""

    GNOMAD = "gnomAD"
    GOOGLE_CLOUD_PUBLIC_DATASETS = "Google Cloud Public Datasets"
    REGISTRY_OF_OPEN_DATA_ON_AWS = "Registry of Open Data on AWS"
    AZURE_OPEN_DATASETS = "Azure Open Datasets"


def get_default_public_resource_source() -> Union[GnomadPublicResourceSource, str]:
    """
    Get the default source for public gnomAD resources.

    :returns: Default resource source
    """
    default_source_from_env = os.getenv("GNOMAD_DEFAULT_PUBLIC_RESOURCE_SOURCE", None)
    if default_source_from_env:
        # Convert to a GnomadPublicResourceSource enum if possible
        try:
            default_source = GnomadPublicResourceSource(default_source_from_env)
            return default_source
        except ValueError:
            return default_source_from_env

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
