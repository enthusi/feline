"""
The `project_path_config` module centralizes path configurations.
It defines paths for the key directories.

Attributes:
    - DATA_PATH_ROOT (str): Root directory path.
    - DATA_PATH_PDF (str): Directory for PDF files.
    - DATA_PATH_LOOKUP (str): Directory for lookup data.
    - DATA_PATH_PROCESSED (str): Directory for processed data.
    - DATA_PATH_RAW (str): Directory for raw data.

Dependencies:
	- os: For file and system operations.
"""

import os

DATA_PATH_ROOT = os.getcwd()
DATA_PATH_PDF = os.path.join(DATA_PATH_ROOT, "data", "pdf_files")
DATA_PATH_LOOKUP = os.path.join(DATA_PATH_ROOT, "data", "lookup")
DATA_PATH_PROCESSED = os.path.join(DATA_PATH_ROOT, "data", "processed")
DATA_PATH_RAW = os.path.join(DATA_PATH_ROOT, "data", "raw")
DATA_PATH_RUNTIME_FILES = os.path.join(DATA_PATH_ROOT, "data", "runtime_files")
