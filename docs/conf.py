# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

sys.path.insert(0, os.path.abspath(".."))


# -- Project information -----------------------------------------------------

project = "bfieldtools"
copyright = "2019, bfieldtools developers"
author = "bfieldtools developers"

# The full version, including alpha/beta/rc tags
release = "0.1"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.


from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.formatting.unsrt import pages, date
from pybtex.style.formatting import toplevel
from pybtex.style.template import (
    field,
    first_of,
    href,
    join,
    names,
    optional,
    optional_field,
    sentence,
    tag,
    together,
    words,
    node,
)
from pybtex.plugin import register_plugin
from pybtex.richtext import Symbol, Text

# bibtex style
class MyStyle(UnsrtStyle):
    def __init__(self):
        super().__init__()
        self.abbreviate_names = True

    # This could be modified
    def get_article_template(self, e):
        volume_and_pages = first_of[
            # volume and pages, with optional issue number
            optional[
                join[field("volume"), optional["(", field("number"), ")"], ":", pages],
            ],
            # pages only
            words["pages", pages],
        ]
        template = toplevel[
            self.format_names("author"),
            self.format_title(e, "title"),
            sentence[tag("em")[field("journal")], optional[volume_and_pages], date],
            self.format_web_refs(e),
            sentence[optional_field("note")],
        ]
        return template


register_plugin("pybtex.style.formatting", "mystyle", MyStyle)

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx_gallery.gen_gallery",
    "sphinx_bootstrap_theme",
    "sphinxcontrib.bibtex",
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# External imports
autodoc_mock_imports = [
    "numpy",
    "scipy",
    "trimesh",
    "mayavi",
    "quadpy",
    "cvxopt",
    "cvxpy",
    "numba",
    "psutil",
]

# Configure sphinx-gallery

scrapers = ("matplotlib",)
try:
    from mayavi import mlab

    # Do not pop up any mayavi windows while running the
    # examples. These are very annoying since they steal the focus.
    mlab.options.offscreen = True
    # hack to initialize the Mayavi Engine
    mlab.test_plot3d()
    mlab.close()
except Exception:
    pass
else:
    scrapers += ("mayavi",)

sphinx_gallery_conf = {
    "examples_dirs": "../examples",  # path to your example scripts
    "gallery_dirs": "auto_examples",  # path where to save gallery generated examples
    "filename_pattern": ".py",  # which examples are executed for plots etc
    "image_scrapers": scrapers,
    "abort_on_example_error": False,
    "download_section_examples": False,
    "show_memory": True,
    "reference_url": {
        # The module you locally document uses None
        "sphinx_gallery": None,
    },
}


# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/devdocs", None),
    "scipy": ("https://scipy.github.io/devdocs", None),
    "matplotlib": ("https://matplotlib.org", None),
    "numba": ("https://numba.pydata.org/numba-doc/latest", None),
    "joblib": ("https://joblib.readthedocs.io/en/latest", None),
    "mayavi": ("http://docs.enthought.com/mayavi/mayavi", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable", None),
}

numpydoc_show_class_members = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = "bootstrap"

html_theme_options = {
    #    'navbar_title': ' ',  # we replace this with an image
    "source_link_position": "nav",  # default
    "bootswatch_theme": "flatly",  # yeti paper lumen
    "navbar_sidebarrel": False,  # Render the next/prev links in navbar?
    "navbar_pagenav": False,
    "navbar_class": "navbar",
    "bootstrap_version": "3",  # default
    "navbar_links": [
        ("Overview", "overview"),
        ("Installation", "installation"),
        ("Literature", "literature"),
        ("Example gallery", "auto_examples/index"),
        ("API Reference", "source/modules"),
    ],
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = "_static/favicon.ico"


# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = "_static/logo_simple.png"
html_logo = "_static/logo_simple.png"

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False
html_copy_source = False
