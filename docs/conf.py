import sphinx_rtd_theme
from docutils.nodes import author
import datetime as dt
import os
import sys

sys.path.insert(0, os.path.abspath('../src/postprocessing'))
sys.path.insert(0, os.path.abspath('../src/preprocessing'))


extensions = [
    'sphinx_rtd_theme',
	'sphinx.ext.autodoc',
	'sphinx.ext.napoleon',  # If you're using Google-style or NumPy-style docstrings
	'sphinx.ext.viewcode',  # Optional: Adds links to highlighted source code
]

project = 'Feline'
copyright = f'{dt.datetime.now():%Y}, Feline Team'
author = 'Martin Wendt, Marvin Henschel, Oskar Soth'
release = '1.0.0'

html_theme = 'sphinx_rtd_theme'
html_show_sourcelink = False


html_context = {
    "display_github": True,
    "github_user": "enthusi",
    "github_repo": "feline",
    "github_version": "main",
}
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'special-members': '__init__',
    'inherited-members': True,
    'show-inheritance': True,
}
