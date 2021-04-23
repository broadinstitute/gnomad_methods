"""Sphinx helpers for documentation."""

import importlib
import symtable
import types

from sphinx.ext.autosummary import Autosummary, get_documenter
from sphinx.pycode import ModuleAnalyzer


class AutoModuleSummary(Autosummary):
    """Autosummary containing all members of a module."""

    required_arguments = 1

    def run(self):
        module_name = str(self.arguments[0])
        module = importlib.import_module(module_name)

        analyzer = ModuleAnalyzer.for_module(module_name)
        attr_docs = analyzer.find_attr_docs()

        with open(module.__file__) as module_file:
            symbols = symtable.symtable(module_file.read(), module.__file__, "exec")

        members = []

        for name in dir(module):
            member = getattr(module, name)

            # Ignore private members
            if name.startswith("_"):
                continue

            # Ignore imported modules
            if isinstance(member, types.ModuleType):
                continue

            # Ignore members imported from other modules
            member_module_name = getattr(member, "__module__", None)
            if member_module_name is None:
                try:
                    if symbols.lookup(name).is_imported():
                        continue
                except KeyError:
                    continue
            else:
                if member_module_name != module_name:
                    continue

            documenter = get_documenter(self.env.app, member, module)

            # Ignore data items that do not have docstrings
            if documenter.objtype == "data" and ("", name) not in attr_docs:
                continue

            members.append(name)

        # Sort members in the order they appear in source code
        tagorder = analyzer.tagorder
        members.sort(key=lambda name: tagorder.get(name, len(tagorder)))

        self.content = [f"{module_name}.{member}" for member in members]

        return super().run()
