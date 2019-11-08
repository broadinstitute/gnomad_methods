from typing import Union, Optional, Tuple, Dict
import hail as hl


# Resource classes
class ResourceBundle:
    def __init__(self, **kwargs):
        self.resources = kwargs
        self.__dict__.update(kwargs)

    def __getitem__(self, item):
        return self.resources[item]

    def __repr__(self):
        return f"ResourceBundle({self.resources.__repr__()})"

    def __str__(self):
        return f"ResourceBundle({self.resources.__str__()})"


class _Resource:
    def __init__(self, path: str, source_path: Optional[str] = None, **kwargs):
        self.source_path = source_path
        self.path = path
        self.__dict__.update(kwargs)

    @classmethod
    def versioned(cls, versions: Dict[str, '_Resource']):
        latest = versions[list(versions)[0]]
        versioned_resource = cls(path=latest.path, source_path=latest.source_path)
        versioned_resource.versions = versions
        return versioned_resource


class TableResource(_Resource):
    def ht(self):
        return hl.read_table(self.path)


class MatrixTableResource(_Resource):
    def mt(self):
        return hl.read_matrix_table(self.path)
