"""Sphinx configuration for the orGUI documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

import orgui

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "orGUI"
copyright = "2026, Timo Fuchs"
author = "Timo Fuchs"
release = orgui.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
]

autodoc_default_options = {
    "members": True,
    "member-order": "bysource",
    "show-inheritance": True,
}
autosummary_generate = True

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_theme_options = {
    "sidebar_hide_name": True,
}
html_logo = "../../orgui/resources/icons/logo_v1.svg"
html_static_path = ["_static", "../../orgui/resources/icons"]
html_css_files = ["custom.css"]
