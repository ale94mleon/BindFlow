# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('./source/'))
sys.path.insert(0, os.path.abspath('./notebooks/'))
from datetime import datetime

# -- Project information -----------------------------------------------------

project = 'MolDrug'
copyright = f"2022-{datetime.now().year}, Alejandro Martínez León"
author = 'Alejandro Martínez León'

# # The full version, including alpha/beta/rc tags
# release = '0.0.1-beta5'


# -- General configuration ---------------------------------------------------

github_doc_root = 'https://github.com/ale94mleon/moldrug/tree/main/docs'
needs_sphinx = '"5.3.0"'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
    "recommonmark",
    "IPython.sphinxext.ipython_console_highlighting",
    "IPython.sphinxext.ipython_directive",
    "myst_nb",
    "sphinx_copybutton",
]


myst_enable_extensions = [
    "colon_fence",
]
nb_execution_allow_errors = False
nb_execution_raise_on_error = True
# nb_execution_timeout = -1
# Do not execute the notebooks, they take too much time
nb_execution_mode = 'off'
myst_heading_anchors = 6

mathjax_path = 'https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML'
autosectionlabel_prefix_document = True
napoleon_google_docstring = True

# copybutton
copybutton_exclude = ".linenos, .gp"

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
pygments_style = "sphinx"
html_theme_options = {
    "repository_url": "https://github.com/ale94mleon/moldrug/",
    "path_to_docs": "docs",
    "use_source_button": True,
    "use_download_button": True,
    "use_repository_button": True,
    "use_issues_button": True,
    "launch_buttons": {"colab_url": "https://colab.research.google.com"},
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/ale94mleon/moldrug/",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        }
    ],
}
html_logo = "source/_static/logo.png"

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [] #['_static']

intersphinx_mapping = {'Python': ('https://docs.python.org/3/', None),
                       'NumPy': ('https://numpy.org/doc/stable/', None),
                       'CReM': ('https://crem.readthedocs.io/en/latest/', None),
                       'RDKit': ('https://www.rdkit.org/docs/', None),
                       'Pandas': ('https://pandas.pydata.org/docs/', None),
                       }