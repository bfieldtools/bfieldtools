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
sys.path.insert(0, os.path.abspath('..'))



# -- Project information -----------------------------------------------------

project = 'bfieldtools'
copyright = '2019, Rasmus Zetter, Antti Mäkinen, Joonas Iivanainen'
author = 'Rasmus Zetter, Antti Mäkinen, Joonas Iivanainen'

# The full version, including alpha/beta/rc tags
release = '0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc',
              'sphinx.ext.coverage',
              'numpydoc',
              'sphinx_gallery.gen_gallery',
              'sphinx_bootstrap_theme']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

#External imports
autodoc_mock_imports = ["numpy", "scipy", "trimesh", "mayavi", "quadpy", "cvxopt", "numba", "psutil"]

#Configure sphinx-gallery

scrapers = ('matplotlib',)
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
    scrapers += ('mayavi',)

sphinx_gallery_conf = {
     'examples_dirs': '../examples',   # path to your example scripts
      'gallery_dirs': 'auto_examples',  # path where to save gallery generated examples
      'filename_pattern':'.py',         # which examples are executed for plots etc
      'image_scrapers': scrapers,
      'abort_on_example_error': False,
      'download_section_examples': False,
      'show_memory': True
}


numpydoc_show_class_members = False

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

html_theme = 'bootstrap'

html_theme_options = {
#    'navbar_title': ' ',  # we replace this with an image
    'source_link_position': "nav",  # default
    'bootswatch_theme': "flatly",  # yeti paper lumen
    'navbar_sidebarrel': False,  # Render the next/prev links in navbar?
    'navbar_pagenav': False,
    'navbar_class': "navbar",
    'bootstrap_version': "3",  # default
    'navbar_links': [
        ("Readme", "readme"),
        ("Literature", "literature"),
        ("Example gallery", "auto_examples/index"),
        ("API Reference", "source/modules")
    ],
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
