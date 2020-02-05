from typing import Optional, Dict, Any, Callable, List
from gnomad_hail.utils.gnomad_functions import logger
import hail as hl
from hail.linalg import BlockMatrix
from abc import ABC, abstractmethod


# Resource classes
class BaseResource(ABC):
    """
    Generic abstract resource class.

    :param path: The resource path
    :param import_sources: Any sources that are required for the import and need to be kept track of (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the resource. `import_func` will be passed the `import_sources` dictionary as kwargs. `import_func` will be passed the `import_sources` dictionary as kwargs.
    :param expected_file_extensions: A list of all expected file extensions. If path doesn't end with one of these, a warning if emitted.
    """

    @abstractmethod
    def __init__(
            self,
            path: Optional[str] = None,
            import_sources: Optional[Dict[str, Any]] = None,
            import_func: Optional[Callable] = None,
            expected_file_extensions: Optional[List[str]] = None
    ):
        if path is None and import_func is None:
            raise ValueError(f"{self.__class__.__name__} requires at least one of path or import_func arguments.")

        self.path = path
        self.import_sources = import_sources
        self.import_func = import_func

        if path is not None and expected_file_extensions and not [ext for ext in expected_file_extensions if path.endswith(ext)]:
            logger.warning(
                "Created the following {} with a path that doesn't ends with {}: {}".format(
                    self.__class__.__name__,
                    " or ".join(expected_file_extensions),
                    self
                )
            )

    def __repr__(self):
        attr_str = [f'path={self.path}']
        if self.import_sources is not None:
            attr_str.append(f'import_sources={self.import_sources}')
        return f'{self.__class__.__name__}({",".join(attr_str)})'

    @abstractmethod
    def import_resource(self):
        """
        Abstract method to import the resource using its import_func and writes it in its path.
        """
        pass


class TableResource(BaseResource):
    """
    A Hail Table resource

    :param path: The Table path (typically ending in .ht)
    :param import_sources: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the Table. `import_func` will be passed the `import_sources` dictionary as kwargs.
    """

    def __init__(
            self,
            path: Optional[str] = None,
            import_sources: Optional[Dict[str, Any]] = None,
            import_func: Optional[Callable[..., hl.Table]] = None
    ):
        super().__init__(
            path=path,
            import_sources=import_sources,
            import_func=import_func,
            expected_file_extensions=['.ht']
        )

    def ht(self, force_import: bool = False) -> hl.Table:
        """
        Read and return the Hail Table resource

        :return: Hail Table resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_sources)
        else:
            return hl.read_table(self.path)

    def import_resource(self, overwrite: bool = False, **kwargs):
        """
        Imports the TableResource using its import_func and writes it in its path.

        :param overwrite: If ``True``, overwrite an existing file at the destination.
        :param kwargs: Any other parameters to be passed to hl.Table.write
        :return: Nothing
        """
        self.import_func(**self.import_sources).write(
            self.path,
            overwrite=overwrite,
            **kwargs
        )


class MatrixTableResource(BaseResource):
    """
    A Hail MatrixTable resource

    :param path: The MatrixTable path (typically ending in .mt)
    :param import_sources: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the MatrixTable. `import_func` will be passed the `import_sources` dictionary as kwargs.
    """

    def __init__(
            self,
            path: Optional[str] = None,
            import_sources: Optional[Dict[str, Any]] = None,
            import_func: Optional[Callable[..., hl.MatrixTable]] = None
    ):
        super().__init__(
            path=path,
            import_sources=import_sources,
            import_func=import_func,
            expected_file_extensions=['.mt']
        )

    def mt(self, force_import: bool = False) -> hl.MatrixTable:
        """
        Read and return the Hail MatrixTable resource

        :return: Hail MatrixTable resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_sources)
        else:
            return hl.read_matrix_table(self.path)

    def import_resource(self, overwrite: bool = False, **kwargs):
        """
        Imports the MatrixTable resource using its import_func and writes it in its path.

        :param overwrite: If set, existing file(s) will be overwritten
        :param kwargs: Any other parameters to be passed to hl.MatrixTable.write
        :return: Nothing
        """
        self.import_func(**self.import_sources).write(self.path, overwrite=overwrite, **kwargs)


class PedigreeResource(BaseResource):
    """
    A pedigree resource

    :param path: The Pedigree path (typically ending in .fam or .ped)
    :param import_sources: Any sources that are required for the import and need to be kept track of and/or passed to the import_func (e.g. .vcf path for an imported VCF)
    :param import_func: A function used to import the MatrixTable. `import_func` will be passed the `import_sources` dictionary as kwargs.
    :param quant_pheno: If ``True``, phenotype is interpreted as quantitative.
    :param delimiter: Field delimiter regex.
    :param missing: The string used to denote missing values. For case-control, 0, -9, and non-numeric are also treated as missing.
    """

    def __init__(
            self,
            path: Optional[str] = None,
            import_sources: Optional[Dict[str, Any]] = None,
            import_func: Optional[Callable[..., hl.MatrixTable]] = None,
            quant_pheno: bool = False,
            delimiter: str = r"\\s+",
            missing: str = 'NA'
    ):
        super().__init__(
            path=path,
            import_sources=import_sources,
            import_func=import_func,
            expected_file_extensions=['.fam', '.ped']
        )

        self.quant_pheno = quant_pheno
        self.delimiter = delimiter
        self.missing = missing

    def ht(self) -> hl.Table:
        """
        Reads the pedigree into a family HT using hl.import_fam().

        :return: Family table
        """
        return hl.import_fam(self.path, quant_pheno=self.quant_pheno, delimiter=self.delimiter, missing=self.missing)

    def pedigree(self) -> hl.Pedigree:
        """
        Reads the pedigree into an hl.Pedigree using hl.Pedigree.read().

        :param delimiter: Delimiter used in the ped file
        :return: pedigree
        """
        return hl.Pedigree.read(self.path, delimiter=self.delimiter)

    def import_resource(self):
        """
        Imports the Pedigree resource using its import_func and writes it in its path.

        :return: Nothing
        """
        self.import_func(**self.import_sources).write(self.path)


class BlockMatrixResource(BaseResource):
    """
    A Hail BlockMatrix resource

    :param path: The BlockMatrix path (typically ending in .bm)
    :param import_sources: Any sources that are required for the import and need to be kept track of and/or passed to the import_func.
    :param import_func: A function used to import the BlockMatrix. `import_func` will be passed the `import_sources` dictionary as kwargs.
    """

    def __init__(
            self,
            path: Optional[str] = None,
            import_sources: Optional[Dict[str, Any]] = None,
            import_func: Optional[Callable[..., BlockMatrix]] = None
    ):
        super().__init__(
            path=path,
            import_sources=import_sources,
            import_func=import_func,
            expected_file_extensions=[".bm"]
        )

    def bm(self) -> BlockMatrix:
        """
        Read and return the Hail MatrixTable resource

        :return: Hail MatrixTable resource
        """
        return BlockMatrix.read(self.path)

    def import_resource(self, overwrite: bool = False, **kwargs):
        """
        Imports the BlockMatrixResource using its import_func and writes it in its path.

        :param overwrite: If ``True``, overwrite an existing file at the destination.
        :param kwargs: Any additional parameters to be passed to BlockMatrix.write
        :return: Nothing
        """
        self.import_func(**self.import_sources).write(self.path, overwrite=False, **kwargs)


class BaseVersionedResource(BaseResource, ABC):
    """
    Abstract class for a versioned resource

    The `path`/`source_path` attributes of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> resource.
    """

    def __init__(self, default_version: str, versions: Dict[str, BaseResource]):
        if type(self) is BaseVersionedResource:
            raise TypeError("Can't instantiate abstract class BaseVersionedResource")

        if default_version not in versions:
            raise KeyError(
                f"default_version {default_version} not found in versions dictionary passed to {self.__class__.__name__}.")

        for version_name, version_resource in versions.items():
            if version_resource.__class__ not in self.__class__.__bases__:
                raise TypeError(f"Cannot create a {self.__class__.__name__} resource with version {version_name} of type {version_resource.__class__.__name__}")

        self.default_version = default_version
        self.versions = versions

        super().__init__(
            path=versions[default_version].path,
            import_sources=versions[default_version].import_sources
        )


    def __repr__(self):
        return f'{self.__class__.__name__}(default_version={self.default_version}, default_resource={self.versions[self.default_version]}, versions={list(self.versions.keys())})'


class VersionedTableResource(BaseVersionedResource, TableResource):
    """
    Versioned Table resource

    The `path`/`source_path` attributes of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this Table resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> TableResource.
    """
    def __init__(self, default_version: str, versions: Dict[str, TableResource]):
        super().__init__(default_version, versions)


class VersionedMatrixTableResource(BaseVersionedResource, MatrixTableResource):
    """
    Versioned MatrixTable resource

    The `path`/`source_path` attributes of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this MatrixTable resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> MatrixTableResource.
    """

    def __init__(self, default_version: str, versions: Dict[str, MatrixTableResource]):
        super().__init__(default_version, versions)


class VersionedPedigreeResource(BaseVersionedResource, PedigreeResource):
    """
    Versioned Pedigree resource

    The `path`/`source_path` attributes of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this Pedigree resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> PedigreeResource.
    """

    def __init__(self, default_version: str, versions: Dict[str, PedigreeResource]):
        super().__init__(default_version, versions)


class VersionedBlockMatrixResource(BaseVersionedResource, BlockMatrixResource):
    """
    Versioned BlockMatrix resource

    The `path`/`source_path` attributes of the versioned resource are those of the default version of the resource.
    In addition, all versions of the resource are stored in the `versions` attribute.

    :param default_version: The default version of this BlockMatrix resource (must to be in the `versions` dict)
    :param versions: A dict of version name -> BlockMatrixResource.
    """

    def __init__(self, default_version: str, versions: Dict[str, BlockMatrixResource]):
        super().__init__(default_version, versions)


class DataException(Exception):
    pass
