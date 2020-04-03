import itertools
import textwrap
from inspect import getmembers
from typing import Dict, Tuple, Optional
from gnomad.resources.resource_utils import BaseResource, BaseVersionedResource
import gnomad.resources.grch37 as grch37
import gnomad.resources.grch38 as grch38
import argparse


# Generate a dictionary of resource available for import for a given genome build
def get_module_importable_resources(
    module, prefix: Optional[str] = None
) -> Dict[str, Tuple[str, BaseResource]]:
    """
    This takes a module that was imported and generates a list of all resources in this module that can be imported (i.e. with a path and import_func).
    The dict produced is as follows:
    keys: {prefix}.{resource_name}.{version} (with prefix only present if `prefix` is set, and `version` only present for versioned resources)
    values: ({resource_name}[ version {version}], resource) with resource_name set to the variable name in the module and the version present for versioned resources.

     The following example will generate a dict with all the resources in gnomad.resources.grch37 that can be imported:

    .. code-block:: python

        import gnomad.resources.grch37 as grch37
        grch37_resources = get_module_importable_resources(grch37, prefix='grch37')
    
    :param module: Input module
    :param prefix: 
    :return: 
    """
    _prefix = f"{prefix}." if prefix else ""
    resources = {}
    for resource_name, resource in getmembers(
        module, lambda x: isinstance(x, BaseResource)
    ):
        if resource.path and resource.import_func:
            arg_name = f"{_prefix}{resource_name}"
            if isinstance(resource, BaseVersionedResource):
                for version in resource.versions:
                    arg_name += f".{version}"
                    resource_name = f"{resource_name} version {version}"
            resources[arg_name] = (resource_name, resource)
    return resources


def get_resources_descriptions(
    resources: Dict[str, Tuple[str, BaseResource]], width: Optional[int] = 100
) -> str:
    """
    Returns a string listing all resources in the input dict along with the path from which they
    are imported and the path at which they are stored.

    :param resources: A dict returned from get_module_importable_resources
    :param width: Maximum width of lines in the returned string
    """
    wrapper = textwrap.TextWrapper(
        width=width, initial_indent=" " * 2, subsequent_indent=" " * 4
    )
    return "\n".join(
        itertools.chain.from_iterable(
            [
                f"{resource_arg}:",
                wrapper.fill(
                    f"import {getattr(resource, 'import_args', {}).get('path', '???')}"
                ),
                wrapper.fill(f"to {resource.path}"),
                "",
            ]
            for resource_arg, (resource_name, resource) in resources.items()
        )
    )


grch37_resources = get_module_importable_resources(grch37, "grch37")
grch38_resources = get_module_importable_resources(grch38, "grch38")
all_resources = {**grch37_resources, **grch38_resources}


def main(args):
    for resource_arg in args.resources:
        resource_name, resource = all_resources[resource_arg]
        print(f"Importing {resource_name}...")
        resource.import_resource(args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        "resources",
        choices=list(all_resources.keys()),
        metavar="resource",
        nargs="+",
        help="Resource to import. Choices are:\n\n"
        + get_resources_descriptions(all_resources),
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    main(parser.parse_args())
