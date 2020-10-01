import typing
from enum import Enum


class GnomadResourceProvider(Enum):
    GNOMAD = "gnomAD"
    GOOGLE_CLOUD_PUBLIC_DATASETS = "Google Cloud Public Datasets"


class GnomadResourceConfiguration:
    """
    Configuration for gnomAD resources.
    """

    __default_resource_provider = GnomadResourceProvider.GOOGLE_CLOUD_PUBLIC_DATASETS

    @property
    def default_resource_provider(self) -> GnomadResourceProvider:
        """
        Get the default provider for resource files.

        This is used to determine resource URLs when `resources_root` is not specified.

        :returns: Provider name
        """
        return self.__default_resource_provider

    @default_resource_provider.setter
    def default_resource_provider(
        self, provider: typing.Union[GnomadResourceProvider, str]
    ) -> None:
        """
        Set the default provider for resource files.

        This is used to determine resource URLs when `resources_root` is not specified.

        :param provider: Provider name
        """
        if not isinstance(provider, GnomadResourceProvider):
            provider = GnomadResourceProvider(provider)

        self.__default_resource_provider = provider


gnomad_resource_configuration = GnomadResourceConfiguration()
