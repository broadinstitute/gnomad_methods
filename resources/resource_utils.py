from typing import Optional, Dict, Any
from gnomad_hail import logger
import hail as hl
from abc import ABC, abstractmethod


# Resource classes
class BaseResource(ABC):
    """
    Generic abstract resource class.
    """

    @abstractmethod
    def __init__(self, path: str, source_path: Optional[str] = None, attributes: Dict[str, Any] = None, expected_file_extension: str = ""):
        """
        Creates a Resource

        :param str path: The resource path
        :param str source_path: The path for the resource source when it was originally imported
        :param dict of str attributes: Additional attributes for the resource
        """

        if expected_file_extension and not path.endswith(expected_file_extension):
            logger.warning(f"Created the following TableResource with a path that doesn't ends with {expected_file_extension}: {self}")

        self.source_path = source_path
        self.path = path
        self.attributes = attributes

    def __repr__(self):
        attr_str = [f'path={self.path}']
        if self.source_path is not None:
            attr_str.append(f'source_path={self.source_path}')
        if self.attributes is not None:
            attr_str.append(f'attributes={self.attributes}')
        return f'{self.__class__.__name__}({",".join(attr_str)})'


class TableResource(BaseResource):
    """
    A Hail Table resource
    """

    def __init__(self, path: str, source_path: Optional[str] = None, attributes: Dict[str, Any] = None):
        super().__init__(
            path=path,
            source_path=source_path,
            attributes=attributes,
            expected_file_extension='.ht'
        )

    def ht(self) -> hl.Table:
        """
        Read and return the Hail Table resource
        :return: Hail Table resource
        :rtype: Table
        """
        return hl.read_table(self.path)


class MatrixTableResource(BaseResource):
    """
    A Hail MatrixTable resource
    """

    def __init__(self, path: str, source_path: Optional[str] = None, attributes: Dict[str, Any] = None):
        super().__init__(
            path=path,
            source_path=source_path,
            attributes=attributes,
            expected_file_extension='.mt'
        )

    def mt(self)-> hl.MatrixTable:
        """
        Read and return the Hail MatrixTable resource
        :return: Hail MatrixTable resource
        :rtype: MatrixTable
        """
        return hl.read_matrix_table(self.path)


class PedigreeResource(BaseResource):
    """
    A pedigree resource
    """

    def __init__(self, path: str, attributes: Dict[str, Any] = None):
        super().__init__(
            path=path,
            source_path=None,
            attributes=attributes,
            expected_file_extension='.ped'
        )

    def ht(self, delimiter=r"\\s+") -> hl.Table:
        """
        Reads the pedigree into a family HT using hl.import_fam().

        :param str delimiter: Delimiter used in the ped file
        :return: Family table
        :rtype: Table
        """
        return hl.import_fam(self.path, delimiter=delimiter)

    def pedigree(self, delimiter=r"\\s+") -> hl.Pedigree:
        """
        Reads the pedigree into an hl.Pedigree using hl.Pedigree.read().

        :param str delimiter: Delimiter used in the ped file
        :return: pedigree
        :rtype: Pedigree
        """
        return hl.Pedigree.read(self.path, delimiter=delimiter)


class BaseVersionedResource(BaseResource, ABC):
    """
    Abstract class for a versioned resource
    """

    def __init__(self, default_version: str, versions: Dict[str, BaseResource]):
        """
        Creates a versioned resource.
        The `path`/`source_path` attributes of the versioned resource are those
        of the default version of the resource.

        In addition, all versions of the resource are stored in the `versions` attribute.

        :param str default_version: The default version of this resource (needs to be in the `versions` dict)
        :param dict of str -> BaseResource versions: A dict of version name -> resource.
        """

        if type(self) is BaseVersionedResource:
            raise TypeError("Can't instantiate abstract class BaseVersionedResource")

        if default_version not in versions:
            raise KeyError(f"default_version {default_version} not found in versions dictionary passed to {self.__class__.__name__}.")

        super().__init__(
            path=versions[default_version].path,
            source_path=versions[default_version].source_path,
            attributes=versions[default_version].attributes
        )
        self.default_version = default_version
        self.versions = versions

    def __repr__(self):
        return f'{self.__class__.__name__}(default_version={self.default_version}, default_resource={self.versions[self.default_version]}, versions={list(self.versions.keys())})'


class VersionedTableResource(TableResource, BaseVersionedResource):
    pass


class VersionedMatrixTableResource(MatrixTableResource, BaseVersionedResource):
    pass
