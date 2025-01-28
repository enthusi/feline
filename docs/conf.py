import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))

extensions = [
    'sphinx_rtd_theme',
	'sphinx.ext.autodoc',
	'sphinx.ext.napoleon',  # If you're using Google-style or NumPy-style docstrings
	'sphinx.ext.viewcode',  # Optional: Adds links to highlighted source code
    'sphinx_copybutton',
]

project = 'Feline'
author = 'Martin Wendt, Marvin Henschel, Oskar Soth'
release = '1.0.0'

html_static_path = ['_static']
html_js_files = [
    'custom.js',  # Add your JavaScript file
]

html_theme = 'sphinx_rtd_theme'
html_show_sourcelink = False
pygments_style = 'sphinx'

autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'private-members': True,
    'special-members': '__init__',
    'inherited-members': True,
    'show-inheritance': True,
}
