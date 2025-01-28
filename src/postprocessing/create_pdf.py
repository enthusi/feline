"""
The `create_pdf` module handles PDF operations.
It provides functionality for merging PDF files and cleaning up the data folder.

Functions:
    - create_pdf_file: Merges all PDF files in the specified folder into one, named with the current timestamp.
    - remove_pdf_files: Deletes all individual PDF files to maintain an organized directory.

Dependencies:
    - pypdf: For merging PDF files.
    - os: For file system operations.
    - datetime: For generating timestamps for the merged PDF.
    - project_path_config: For accessing project-specific path configurations.
"""

import datetime
import os
import pypdf

import src.project_path_config as project_path_config


def create_pdf_file() -> None:
    """
    This function creates a PDF file by merging all single PDF files
    in the data folder.

    Args:
        None

    Returns:
        None
    """
    merger = pypdf.PdfWriter()
    files_in_directory = os.listdir(project_path_config.DATA_PATH_PDF)
    files_in_directory.sort(reverse=True)
    pdf_files = []

    for pdf_file in files_in_directory:
        if "fig" in pdf_file and pdf_file.endswith(".pdf"):
            pdf_files.append(os.path.join(project_path_config.DATA_PATH_PDF,
                                          pdf_file))

    for pdf_file in pdf_files:
        merger.append(pdf_file)
    timestamp = datetime.datetime.now().strftime("%Y_%m_%d_%H:%M:%S")
    file_name = f"result_{timestamp}.pdf"
    merger.write(os.path.join(project_path_config.DATA_PATH_PDF, file_name))
    merger.close()
    print("\033[92m All single PDF files successfully merged! \033[0m")


def remove_pdf_files() -> None:
    """
    This function deletes all single PDF files in the data folder.

    Args:
        None

    Returns:
        None
    """
    pdf_file_list = os.listdir(project_path_config.DATA_PATH_PDF)
    for pdf_file in range(len(pdf_file_list)):
        if "fig" in pdf_file_list[pdf_file]:
            os.remove(os.path.join(project_path_config.DATA_PATH_PDF,
                                   pdf_file_list[pdf_file]))


# print("\033[91m All single PDF files successfully deleted! \033[0m")


if __name__ == "__main__":
    create_pdf_file()
    remove_pdf_files()
