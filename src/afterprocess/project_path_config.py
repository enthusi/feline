import os

current_level_directory = os.getcwd()
parent_level_directory = os.path.abspath(os.path.join(current_level_directory, os.pardir))
top_level_directory = os.path.abspath(os.path.join(parent_level_directory, os.pardir))

PROJECT_ROOT = top_level_directory

DATA_PATH_RAW = "../../data/raw/"
DATA_PATH_PROCESSED = "../../data/processed/"
