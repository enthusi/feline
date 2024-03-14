import os
import pypdf
import project_path_config


def create_pdf_file():
	merger = pypdf.PdfWriter()
	files_in_directory = os.listdir(project_path_config.DATA_PATH_PDF)
	files_in_directory.sort(reverse=True)
	pdf_files = []

	for pdf_file in files_in_directory:
		if pdf_file.endswith(".pdf"):
			pdf_files.append(os.path.join(project_path_config.DATA_PATH_PDF, pdf_file))

	for pdf_file in pdf_files:
		merger.append(pdf_file)

	merger.write(os.path.join(project_path_config.DATA_PATH_PDF, "all_images.pdf"))
	merger.close()
	print("\033[92m All single PDF files successfully merged! \033[0m")


def remove_pdf_files():
	pdf_file_list = os.listdir(project_path_config.DATA_PATH_PDF)
	for pdf_file in range(len(pdf_file_list)):
		if pdf_file_list[pdf_file].endswith(".pdf") and not pdf_file_list[pdf_file] == "all_images.pdf":
			os.remove(os.path.join(project_path_config.DATA_PATH_PDF, pdf_file_list[pdf_file]))
	# print("\033[91m All single PDF files successfully deleted! \033[0m")


create_pdf_file()
remove_pdf_files()
