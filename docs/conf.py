import sys
from pathlib import Path

from directives import AutoModuleSummary

# Add gnomad_hail to import path.
sys.path.insert(0, str(Path(__file__).parent.parent))


project = "gnomad_hail"
version = release = "master"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
]

master_doc = "index"

html_theme = "sphinx_rtd_theme"

html_theme_options = {
    "display_version": False,
}

html_static_path = ["_static"]

html_css_files = ["theme_overrides.css"]

html_show_sphinx = False

html_show_copyright = False

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "member-order": "bysource",
}


def setup(app):
    app.add_directive("gnomadhail_automodulesummary", AutoModuleSummary)
