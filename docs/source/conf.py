# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys
sys.path.insert(0, os.path.abspath('../../src/glidesign'))

import sphinx_rtd_theme


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GliDesign'
copyright = '2025, Ties Leniger, Niels van Nieuwland'
author = 'Ties Leniger, Niels van Nieuwland'
release = '0.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary']
autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

# -- Options for latex output
latex_documents = [
    ('index', 'glidesign.tex', 'GliDesign Documentation', 'Ties Leniger, Niels van Nieuwland', 'manual'),
]

latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    'preamble': r'''
        \usepackage[utf8]{inputenc}       % Encoding
        \usepackage{amsmath}              % Math packages
        \usepackage{amsfonts}             % Fonts for math symbols
        \usepackage{amssymb}              % Extra symbols
        \usepackage{graphicx}             % For including images
        \usepackage{hyperref}             % Links within document and external links
        \usepackage{geometry}             % Page layout and margins
        \usepackage{fancyhdr}             % Custom headers/footers
        \usepackage{longtable}            % For tables spanning multiple pages
        \usepackage{listings}             % For including code blocks
        \usepackage{color}                % For text coloring
        \usepackage{tikz}                 % For drawing graphics (if needed)
        \usepackage{subcaption}           % For subfigures
        \usepackage{multirow}             % For multirow in tables
        \usepackage{pdfpages}             % To include PDFs in the document (if necessary)
        \usepackage{parskip}              % To adjust paragraph skipping
        \usepackage{setspace}             % Line spacing control
        \usepackage{tocbibind}            % To include bibliography in the TOC
        \usepackage{hyperref}             % Hyperlinks in PDF
        \hypersetup{colorlinks=true, linkcolor=blue, filecolor=magenta, urlcolor=blue, pdftitle={GliDesign Documentation}, pdfpagemode=FullScreen, pdfauthor={Ties Leniger, Niels van Nieuwland}}
    ''',
    'figure_align': 'htbp',
}