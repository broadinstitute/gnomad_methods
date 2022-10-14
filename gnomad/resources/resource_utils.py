# noqa: D100

import logging
from abc import ABC, abstractmethod
from functools import reduce, wraps
from typing import Any, Callable, Dict, Iterable, List, Optional

import hail as hl
from hail.linalg import BlockMatrix

from gnomad.resources.config import (
    GnomadPublicResourceSource,
    gnomad_public_resource_configuration,
)

logger = logging.getLogger("gnomad.resources")


GNOMAD_PUBLIC_BUCKETS = ("gnomad-public", "gnomad-public-requester-pays")
"""
Public buckets used to stage gnomAD data.

`gnomad-public` is a legacy bucket and contains one readme text file.

The gnomAD Production Team writes output data to `gnomad-public-requester-pays`, and all data in this bucket
syncs to the public bucket `gcp-public-data--gnomad`.
"""

# Resource classes


class BaseResource(ABC):
    """
    Generic abstract resource class.

    :param path: The resource path
    :param import_args: Any sources that are required for the import and need to be kept track of (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the resource. `import_func` will be passed the `import_args` dictionary as kwargs.
    """

    expected_file_extensions: List[str] = []
    """Expected file extensions for this resource type. If path doesn't end with one of these, a warning is logged."""

    def __init__(
        self,
        path: Optional[str] = None,
        import_args: Optional[Dict[str, Any]] = None,
        import_func: Optional[Callable] = None,
    ):
        if path is None and import_func is None:
            raise ValueError(
                f"{self.__class__.__name__} requires at least one of path or"
                " import_func arguments."
            )

        self.path = path
        self.import_args = import_args
        self.import_func = import_func

        if (
            path is not None
            and self.expected_file_extensions
            and not any(path.endswith(ext) for ext in self.expected_file_extensions)
        ):
            logger.warning(
                "Created the following %s with a path that doesn't end with %s: %s",
                self.__class__.__name__,
                " or ".join(self.expected_file_extensions),
                self,
            )

    def __repr__(self):
        attr_str = [f"path={self._path}"]
        if self.import_args is not None:
            attr_str.append(f"import_args={self.import_args}")
        return f'{self.__class__.__name__}({",".join(attr_str)})'

    def _get_path(self):
        return self._path

    def _set_path(self, path):
        self._path = path  # pylint: disable=attribute-defined-outside-init

    # Defining path property this way instead of using a decorator allows _get_path and _set_path
    # to be overridden in subclasses without having to reconfigure the property.
    path = property(
        fget=lambda self: self._get_path(),
        fset=lambda self, path: self._set_path(path),
    )

    @abstractmethod
    def import_resource(self, overwrite: bool = True, **kwargs) -> None:
        """
        Abstract method to import the resource using its import_func and writes it in its path.

        :param overwrite: If ``True``, overwrite an existing file at the destination.
        :param kwargs: Any other parameters to be passed to the underlying hail write function (acceptable parameters depend on specific resource types)
        """


class TableResource(BaseResource):
    """
    A Hail Table resource.

    :param path: The Table path (typically ending in .ht)
    :param import_args: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the Table. `import_func` will be passed the `import_args` dictionary as kwargs.
    """

    expected_file_extensions: List[str] = [".ht"]

    def ht(self, force_import: bool = False) -> hl.Table:
        """
        Read and return the Hail Table resource.

        :return: Hail Table resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_args)
        else:
            return hl.read_table(self.path)

    def import_resource(self, overwrite: bool = True, **kwargs) -> None:
        """
        Import the TableResource using its import_func and writes it in its path.

        :param overwrite: If ``True``, overwrite an existing file at the destination.
        :param kwargs: Any other parameters to be passed to hl.Table.write
        :return: Nothing
        """
        self.import_func(**self.import_args).write(
            self.path, overwrite=overwrite, **kwargs
        )


class MatrixTableResource(BaseResource):
    """
    A Hail MatrixTable resource.

    :param path: The MatrixTable path (typically ending in .mt)
    :param import_args: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the MatrixTable. `import_func` will be passed the `import_args` dictionary as kwargs.
    """

    expected_file_extensions: List[str] = [".mt"]

    def mt(self, force_import: bool = False) -> hl.MatrixTable:
        """
        Read and return the Hail MatrixTable resource.

        :return: Hail MatrixTable resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_args)
        else:
            return hl.read_matrix_table(self.path)

    def import_resource(self, overwrite: bool = True, **kwargs) -> None:
        """
        Import the MatrixTable resource using its import_func and writes it in its path.

        :param overwrite: If set, existing file(s) will be overwritten
        :param kwargs: Any other parameters to be passed to hl.MatrixTable.write
        :return: Nothing
        """
        self.import_func(**self.import_args).write(
            self.path, overwrite=overwrite, **kwargs
        )


class VariantDatasetResource(BaseResource):
    """
    A Hail VariantDataset resource.

    :param path: The VariantDataset path (typically ending in .vds)
    :param import_args: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the VariantDataset. `import_func` will be passed the `import_args` dictionary as kwargs.
    """

    expected_file_extensions: List[str] = [".vds"]

    def vds(self, force_import: bool = False) -> hl.vds.VariantDataset:
        """
        Read and return the Hail VariantDataset resource.

        :return: Hail VariantDataset resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_args)
        else:
            return hl.vds.read_vds(self.path)

    def import_resource(self, overwrite: bool = True, **kwargs) -> None:
        """
        Import the VariantDataset resource using its import_func and writes it in its path.

        :param overwrite: If set, existing file(s) will be overwritten
        :param kwargs: Any other parameters to be passed to hl.vds.VariantDataset.write
        :return: Nothing
        """
        self.import_func(**self.import_args).write(
            self.path, overwrite=overwrite, **kwargs
        )


class PedigreeResource(BaseResource):
    """
    A pedigree resource.

    :param path: The Pedigree path (typically ending in .fam or .ped)
    :param import_args: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the Pedigree. `import_func` will be passed the `import_args` dictionary as kwargs.
    :param quant_pheno: If ``True``, phenotype is interpreted as quantitative.
    :param delimiter: Field delimiter regex.
    :param missing: The string used to denote missing values. For case-control, 0, -9, and non-numeric are also treated as missing.
    """

    expected_file_extensions: List[str] = [".fam", ".ped"]

    def __init__(
        self,
        path: Optional[str] = None,
        import_args: Optional[Dict[str, Any]] = None,
        import_func: Optional[Callable[..., hl.Pedigree]] = None,
        quant_pheno: bool = False,
        delimiter: str = r"\\s+",
        missing: str = "NA",
    ):
        super().__init__(
            path=path,
            import_args=import_args,
            import_func=import_func,
        )

        self.quant_pheno = quant_pheno
        self.delimiter = delimiter
        self.missing = missing

    def ht(self) -> hl.Table:
        """
        Read the pedigree into a family HT using hl.import_fam().

        :return: Family table
        """
        return hl.import_fam(
            self.path,
            quant_pheno=self.quant_pheno,
            delimiter=self.delimiter,
            missing=self.missing,
        )

    def pedigree(self) -> hl.Pedigree:
        """
        Read the pedigree into an hl.Pedigree using hl.Pedigree.read().

        :param delimiter: Delimiter used in the ped file
        :return: pedigree
        """
        return hl.Pedigree.read(self.path, delimiter=self.delimiter)

    def import_resource(self, overwrite: bool = True, **kwargs) -> None:
        """
        Import the Pedigree resource using its import_func and writes it in its path.

        :param overwrite: If set, existing file(s) will be overwritten. IMPORTANT: Currently there is no implementation of this method when `overwrite` is set the `False`
        :param kwargs: Any other parameters to be passed to hl.Pedigree.write
        :return: Nothing
        """
        if not overwrite:
            raise NotImplementedError

        self.import_func(**self.import_args).write(self.path)


class BlockMatrixResource(BaseResource):
    """
    A Hail BlockMatrix resource.

    :param path: The BlockMatrix path (typically ending in .bm)
    :param import_args: Any sources that are required for the import and need to be kept track of and/or passed to the import_func.
    :param import_func: A function used to import the BlockMatrix. `import_func` will be passed the `import_args` dictionary as kwargs.
    """

    expected_file_extensions: List[str] = [".bm"]

    def bm(self) -> BlockMatrix:
        """
        Read and return the Hail MatrixTable resource.

        :return: Hail MatrixTable resource
        """
        return BlockMatrix.read(self.path)

    def import_resource(self, overwrite: bool = True, **kwargs) -> None:
        """
        Import the BlockMatrixResource using its import_func and writes it in its path.

        :param overwrite: If ``True``, overwrite an existing file at the destination.
        :param kwargs: Any additional parameters to be passed to BlockMatrix.write
        :return: Nothing
        """
        self.import_func(**self.import_args).write(
            self.path, overwrite=overwrite, **kwargs
        )


class BaseVersionedResource:
    """
    Class for a versioned resource.

    The attributes and methods of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this resource (must be in the `versions` dict)
    :param versions: A dict of version name -> resource.
    """

    resource_class = BaseResource

    __slots__ = {"default_version", "versions"}

    def __init__(self, default_version: str, versions: Dict[str, BaseResource]):
        default_resource = versions[default_version]

        for version_resource in versions.values():
            if not isinstance(version_resource, self.resource_class):
                raise TypeError(
                    f"{self.__class__.__name__} requires all versions to be of type"
                    f" {self.resource_class.__name__}"
                )

            if version_resource.__class__ is not default_resource.__class__:
                raise TypeError(
                    f"{self.__class__.__name__} requires all versions to be of the same"
                    " type"
                )

        self.default_version = default_version
        self.versions = versions

    def __repr__(self):
        return (
            "{cls}(default_version={default_version}, versions={{{versions}}})".format(
                cls=self.__class__.__name__,
                default_version=self.default_version,
                versions=", ".join(
                    f'"{k}": {repr(v)}' for k, v in self.versions.items()
                ),
            )
        )

    def __getattr__(self, name):
        # If __getattr__ is called for 'default_version', 'version', etc. then
        # something has gone wrong.
        if name in self.__slots__:
            raise ValueError("VersionedResource has not been initialized")

        return getattr(self.versions[self.default_version], name)


class VersionedTableResource(BaseVersionedResource):
    """
    Versioned Table resource.

    The attributes (path, import_args and import_func) of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this Table resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> TableResource.
    """

    resource_class = TableResource

    def __init__(self, default_version: str, versions: Dict[str, TableResource]):
        super().__init__(default_version, versions)


class VersionedMatrixTableResource(BaseVersionedResource):
    """
    Versioned MatrixTable resource.

    The attributes (path, import_args and import_func) of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this MatrixTable resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> MatrixTableResource.
    """

    resource_class = MatrixTableResource

    def __init__(self, default_version: str, versions: Dict[str, MatrixTableResource]):
        super().__init__(default_version, versions)


class VersionedVariantDatasetResource(BaseVersionedResource):
    """
    Versioned VariantDataset resource.

    The attributes (path, import_args and import_func) of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.
    :param default_version: The default version of this VariantDataset resource (must to be in the `versions` dict)

    :param versions: A dict of version name -> VariantDatasetResource.
    """

    resource_class = VariantDatasetResource

    def __init__(
        self, default_version: str, versions: Dict[str, VariantDatasetResource]
    ):
        super().__init__(default_version, versions)


class VersionedPedigreeResource(BaseVersionedResource, PedigreeResource):
    """
    Versioned Pedigree resource.

    The attributes (path, import_args and import_func) of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this Pedigree resource (must be in the `versions` dict)
    :param versions: A dict of version name -> PedigreeResource.
    """

    resource_class = PedigreeResource

    def __init__(self, default_version: str, versions: Dict[str, PedigreeResource]):
        super().__init__(default_version, versions)


class VersionedBlockMatrixResource(BaseVersionedResource, BlockMatrixResource):
    """
    Versioned BlockMatrix resource.

    The attributes (path, import_args and import_func) of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this BlockMatrix resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> BlockMatrixResource.
    """

    resource_class = BlockMatrixResource

    def __init__(self, default_version: str, versions: Dict[str, BlockMatrixResource]):
        super().__init__(default_version, versions)


class ResourceNotAvailable(Exception):
    """Exception raised if a resource is not available from the selected source."""


class GnomadPublicResource(BaseResource, ABC):
    """Base class for the gnomAD project's public resources."""

    def __init_subclass__(cls, *, read_resource_methods: Iterable[str] = []) -> None:
        super().__init_subclass__()

        # Some resources may not be available from all sources due to delays in syncing, etc.
        # This wraps all methods that read the resource and adds a check for if the resource
        # is available from the selected source. If the resource is not available, this
        # throws a more helpful error than if the read were simply allowed to fail.
        def _wrap_read_resource_method(method_name):
            original_method = getattr(cls, method_name)

            @wraps(original_method)
            def read_resource(self, *args, **kwargs):
                # If one of the known sources is selected, check if the resource is available.
                # For custom sources, skip the check and attempt to read the resource.
                resource_source = gnomad_public_resource_configuration.source
                if not self.is_resource_available():
                    if resource_source == GnomadPublicResourceSource.GNOMAD:
                        message = (
                            "This resource is not currently available from the gnomAD"
                            " project public buckets."
                        )
                    elif isinstance(resource_source, GnomadPublicResourceSource):
                        message = (
                            "This resource is not currently available from"
                            f" {resource_source.value}."
                        )
                    else:
                        message = (
                            "This resource is not currently available from"
                            f" {resource_source}."
                        )

                    raise ResourceNotAvailable(
                        f"{message}\n\nTo load resources from a different source (for"
                        " example, Google Cloud Public Datasets) instead, use:\n\n>>>"
                        " from gnomad.resources.config import"
                        " gnomad_public_resource_configuration,"
                        " GnomadPublicResourceSource\n>>>"
                        " gnomad_public_resource_configuration.source ="
                        " GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS\n\nTo"
                        " get all available sources for gnomAD resources, use:\n\n>>>"
                        " from gnomad.resources.config import"
                        " GnomadPublicResourceSource\n>>>"
                        " list(GnomadPublicResourceSource)"
                    )

                return original_method(self, *args, **kwargs)

            setattr(cls, method_name, read_resource)

        for method_name in read_resource_methods:
            _wrap_read_resource_method(method_name)

    def _get_path(self) -> str:
        resource_source = gnomad_public_resource_configuration.source
        if resource_source == GnomadPublicResourceSource.GNOMAD:
            return self._path

        relative_path = reduce(
            lambda path, bucket: path[5 + len(bucket) :]
            if path.startswith(f"gs://{bucket}/")
            else path,
            GNOMAD_PUBLIC_BUCKETS,
            self._path,
        )

        if resource_source == GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS:
            return f"gs://gcp-public-data--gnomad{relative_path}"

        if resource_source == GnomadPublicResourceSource.REGISTRY_OF_OPEN_DATA_ON_AWS:
            return f"s3a://gnomad-public-us-east-1{relative_path}"

        if resource_source == GnomadPublicResourceSource.AZURE_OPEN_DATASETS:
            return f"wasbs://dataset@datasetgnomad.blob.core.windows.net{relative_path}"

        return (
            f"{resource_source.rstrip('/')}{relative_path}"  # pylint: disable=no-member
        )

    def _set_path(self, path):
        if not any(
            path.startswith(f"gs://{bucket}/") for bucket in GNOMAD_PUBLIC_BUCKETS
        ):
            raise ValueError(
                "GnomadPublicResource requires a path to a file in one of the public"
                f" gnomAD buckets ({', '.join(GNOMAD_PUBLIC_BUCKETS)})"
            )

        return super()._set_path(path)

    def is_resource_available(self) -> bool:
        """
        Check if this resource is available from the selected source.

        :return: True if the resource is available.
        """
        path = self.path

        # Hail Tables, MatrixTables, and BlockMatrices are directories.
        # For those, check for the existence of the _SUCCESS object.
        path_to_test = (
            f"{path}/_SUCCESS"
            if any(path.endswith(ext) for ext in (".ht", ".mt", ".bm"))
            else path
        )

        return hl.current_backend().fs.exists(path_to_test)


class GnomadPublicTableResource(
    TableResource, GnomadPublicResource, read_resource_methods=("ht",)
):
    """Resource class for a public Hail Table published by the gnomAD project."""


class GnomadPublicMatrixTableResource(
    MatrixTableResource, GnomadPublicResource, read_resource_methods=("mt",)
):
    """Resource class for a public Hail MatrixTable published by the gnomAD project."""


class GnomadPublicPedigreeResource(
    PedigreeResource, GnomadPublicResource, read_resource_methods=("ht", "pedigree")
):
    """Resource class for a public pedigree published by the gnomAD project."""


class GnomadPublicBlockMatrixResource(
    BlockMatrixResource, GnomadPublicResource, read_resource_methods=("bm",)
):
    """Resource class for a public Hail BlockMatrix published by the gnomAD project."""


class DataException(Exception):  # noqa: D101
    pass


NO_CHR_TO_CHR_CONTIG_RECODING = {
    "1": "chr1",
    "2": "chr2",
    "3": "chr3",
    "4": "chr4",
    "5": "chr5",
    "6": "chr6",
    "7": "chr7",
    "8": "chr8",
    "9": "chr9",
    "10": "chr10",
    "11": "chr11",
    "12": "chr12",
    "13": "chr13",
    "14": "chr14",
    "15": "chr15",
    "16": "chr16",
    "17": "chr17",
    "18": "chr18",
    "19": "chr19",
    "20": "chr20",
    "21": "chr21",
    "22": "chr22",
    "X": "chrX",
    "Y": "chrY",
    "MT": "chrM",
}

DBSNP_B154_CHR_CONTIG_RECODING = {
    "NC_000001.11": "chr1",
    "NC_000002.12": "chr2",
    "NC_000003.12": "chr3",
    "NC_000004.12": "chr4",
    "NC_000005.10": "chr5",
    "NC_000006.12": "chr6",
    "NC_000007.14": "chr7",
    "NC_000008.11": "chr8",
    "NC_000009.12": "chr9",
    "NC_000010.11": "chr10",
    "NC_000011.10": "chr11",
    "NC_000012.12": "chr12",
    "NC_000013.11": "chr13",
    "NC_000014.9": "chr14",
    "NC_000015.10": "chr15",
    "NC_000016.10": "chr16",
    "NC_000017.11": "chr17",
    "NC_000018.10": "chr18",
    "NC_000019.10": "chr19",
    "NC_000020.11": "chr20",
    "NC_000021.9": "chr21",
    "NC_000022.11": "chr22",
    "NC_000023.11": "chrX",
    "NC_000024.10": "chrY",
}


def import_sites_vcf(**kwargs) -> hl.Table:
    """Import site-level data from a VCF into a Hail Table."""
    return hl.import_vcf(**kwargs).rows()
