"""
This module, `project_path_config`, centralizes the path configurations for an
astronomical data analysis project. It defines the root directory and specific
subdirectories for storing PDF files, lookup data, processed data,
and raw data.

Attributes:
    DATA_PATH_ROOT (str): The absolute path to the project's root directory.
    DATA_PATH_PDF (str): The path to the directory for storing PDF files.
    DATA_PATH_LOOKUP (str): The path to the directory for lookup data.
    DATA_PATH_PROCESSED (str): The path to the directory for processed data.
    DATA_PATH_RAW (str): The path to the directory for raw data.

This configuration module ensures that paths are consistently managed across
different parts of the project, facilitating access to data files and
directories. It uses the `os` module to construct absolute paths, making it
easier to reference data across the project's codebase.

Note:
This module is essential for maintaining an organized file structure and
should be imported wherever file paths are needed within the project.
"""

import os

DATA_PATH_ROOT = os.path.abspath("../../")
DATA_PATH_PDF = os.path.join(DATA_PATH_ROOT, "data", "pdf_files")
DATA_PATH_LOOKUP = os.path.join(DATA_PATH_ROOT, "data", "lookup")
DATA_PATH_PROCESSED = os.path.join(DATA_PATH_ROOT, "data", "processed")
DATA_PATH_RAW = os.path.join(DATA_PATH_ROOT, "data", "raw")
DATA_PATH_RUNTIME_FILES = os.path.join(DATA_PATH_ROOT, "data", "runtime_files")
