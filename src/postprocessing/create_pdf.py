"""
This module, `create_pdf`, is designed for handling PDF operations within the
astronomical data analysis project. It provides functionalities to merge single
PDF files into a comprehensive document and to clean up the data folder by
removing these single PDF files after merging. This is particularly useful for
aggregating individual plot files into a single report document and maintaining
a clean workspace.

Functions:
- create_pdf_file: Merges all single PDF files found in the specified data
  folder into one PDF file, named with the current timestamp to
  ensure uniqueness.
- remove_pdf_files: Deletes all single PDF files in the data folder to free up
  space and keep the directory organized.

The module leverages the `pypdf` library for PDF manipulation,
including merging PDF documents. It also uses configurations from
`project_path_config` to locate the data folder where PDF files are stored.
This module is part of a larger project focused on the analysis of astronomical
data, facilitating the documentation and reporting aspect of the analysis.

Dependencies:
- pypdf: For PDF file manipulation, specifically merging PDF files.
- os: For interacting with the file system.
- datetime: For generating timestamps used in naming the merged PDF file.
- project_path_config: For accessing project-specific path configurations.

Note:
This module assumes that the PDF files to be merged follow a specific naming
convention (containing 'fig' in their names) and are located in a predefined
data folder as specified in `project_path_config`.
"""

import os
import pypdf
import project_path_config
import datetime


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
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
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
