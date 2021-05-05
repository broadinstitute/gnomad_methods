"""Configuration for loading resources."""

from enum import Enum
from typing import Union


class GnomadPublicResourceSource(Enum):
    """Sources for public gnomAD resources."""

    GNOMAD = "gnomAD"
    GOOGLE_CLOUD_PUBLIC_DATASETS = "Google Cloud Public Datasets"


DEFAULT_GNOMAD_PUBLIC_RESOURCE_SOURCE = GnomadPublicResourceSource.GNOMAD


class _GnomadPublicResourceConfiguration:
    """Configuration for public gnomAD resources."""

    __source: Union[
        GnomadPublicResourceSource, str
    ] = DEFAULT_GNOMAD_PUBLIC_RESOURCE_SOURCE

    @property
    def source(self) -> Union[GnomadPublicResourceSource, str]:
        """
        Get the source for public gnomAD resource files.

        This is used to determine which URLs gnomAD resources will be loaded from.

        :returns: Source name or path to root of resources directory
        """
        return self.__source

    @source.setter
    def source(self, source: Union[GnomadPublicResourceSource, str]) -> None:
        """
        Set the default source for resource files.

        This is used to determine which URLs gnomAD resources will be loaded from.

        :param source: Source name or path to root of resources directory
        """
        self.__source = source


gnomad_public_resource_configuration = _GnomadPublicResourceConfiguration()
