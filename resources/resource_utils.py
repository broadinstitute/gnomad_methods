from typing import Optional, Dict
import hail as hl


# Resource classes
class _Resource:
    """
    Generic abstract resource class.
    """
    def __init__(self, path: str, source_path: Optional[str] = None, **kwargs):
        if type(self) is _Resource:
            raise TypeError("Can't instantiate abstract class _Resource")

        self.source_path = source_path
        self.path = path
        self.__dict__.update(kwargs)

    @classmethod
    def versioned(cls, versions: Dict[str, '_Resource']) -> '_Resource':
        """
        Creates a versioned resource.
        The `path`/`source_path` attributes of the versioned resource are those
        of default version of the resource (first in the dict keys).

        In addition, all versions of the resource are stored in the `versions` attribute.

        :param dict of str -> _Resource versions: A dict of version name -> resource. The default version should be first in the dict keys.
        :return: A Resource with versions
        :rtype: _Resource
        """
        latest = versions[list(versions)[0]]
        versioned_resource = cls(path=latest.path, source_path=latest.source_path)
        versioned_resource.versions = versions
        return versioned_resource


class TableResource(_Resource):
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


class MatrixTableResource(_Resource):
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
