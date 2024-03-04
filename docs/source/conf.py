# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

import toml

topdir = "../../"

sys.path.insert(0, os.path.abspath(topdir))

pyproject_settings = toml.load("../../pyproject.toml")
project = pyproject_settings["tool"]["poetry"]["name"]

pyproject_settings = toml.load("../../pyproject.toml")
project = pyproject_settings["tool"]["poetry"]["name"]
author = " ,".join(pyproject_settings["tool"]["poetry"]["authors"])
copyright = "CERN"
release = pyproject_settings["tool"]["poetry"]["version"]
version = release


# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.autodoc',
              "sphinx.ext.autosummary",
               "sphinx.ext.viewcode",
               "sphinx.ext.napoleon",
              'sphinx.ext.duration',
              'sphinx.ext.doctest',
              'sphinx_gallery.gen_gallery',
              "sphinxcontrib.apidoc",
              'myst_parser']

apidoc_module_dir = topdir + 'cauchy'
apidoc_output_dir = "gen"
apidoc_separate_modules = True

templates_path = ['_templates']
exclude_patterns = []

source_suffix = ['.rst', '.md']

sphinx_gallery_conf = {'examples_dirs': '../../examples',   # path to your example scripts
     	'gallery_dirs': 'auto_examples',                    # path to where to save gallery generated output
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output


html_theme = 'pydata_sphinx_theme'
# html_static_path = ['_static']
