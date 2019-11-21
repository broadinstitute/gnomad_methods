from typing import Optional, Dict, Any
import hail as hl


# Resource classes
class BaseResource:
    """
    Generic abstract resource class.
    """

    versions = None
    "Versions of this resource"


    def __init__(self, path: str, source_path: Optional[str] = None, attributes: Dict[str, Any] = None):
        """
        Creates a Resource

        :param str path: The resource path
        :param str source_path: The path for the resource source when it was originally imported
        :param dict of str attributes: Additional attributes for the resource
        """
        if type(self) is BaseResource:
            raise TypeError("Can't instantiate abstract class _Resource")

        self.source_path = source_path
        self.path = path
        self.attributes = attributes

    def __repr__(self):
        attr_str = [fr'path={self.path}']
        if self.source_path is not None:
            attr_str.append(fr'source_path={self.source_path}')
        if self.versions is not None:
            attr_str.append(fr'versions={list(self.versions.keys())}')
        if self.attributes is not None:
            attr_str.append(fr'attributes={self.attributes}')
        return f'{self.__class__.__name__}({",".join(attr_str)})'

    @classmethod
    def versioned(cls, default_version: str, versions: Dict[str, 'BaseResource']) -> 'BaseResource':
        """
        Creates a versioned resource.
        The `path`/`source_path` attributes of the versioned resource are those
        of default version of the resource.

        In addition, all versions of the resource are stored in the `versions` attribute.

        :param str default_version: The default version of this resource (needs to be in the `versions` dict)
        :param dict of str -> BaseResource versions: A dict of version name -> resource.
        :return: A Resource with versions
        :rtype: BaseResource
        """
        if default_version not in versions:
            raise KeyError(f"default_version {default_version} not found in versions dictionary passed to {cls.__name__}.versioned()." )

        default_resource = versions[default_version]
        if default_resource.versions is not None:
            if default_resource.versions == versions:
                return default_resource
            else:
                raise ValueError(f"The default resource specified ({default_version}) for creating a versioned resource has a versions attribute that contains different versions than those passed to {cls.__name__}.versioned()")

        default_resource.versions = versions
        return default_resource


class TableResource(BaseResource):
    """
    A Hail Table resource
    """
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
    def mt(self)-> hl.MatrixTable:
        """
        Read and return the Hail MatrixTable resource
        :return: Hail MatrixTable resource
        :rtype: MatrixTable
        """
        return hl.read_matrix_table(self.path)
