import os
import struct
import sys
import mpdaf
import project_path_config


def load_data(file, ext):
    """
    Load data from a file.

    Parameters:
    file (str): The file path.
    ext (int): The extension number.

    Returns:
    mpdaf.obj.Cube: The sum of the cube data along the specified axis.
    """
    return mpdaf.obj.Cube(file, ext=ext).sum(axis=0)


def create_masking_plot():
    """
    Create a masking plot and write it to a file.

    The function takes no parameters. It uses the filename from the command line arguments,
    loads the data and statistic from the file, creates a signal-to-noise image, and writes it to 'image00.fits'.
    """
    filename_first_argument = sys.argv[1]
    file = os.path.join(project_path_config.DATA_PATH_RAW, filename_first_argument)

    data = load_data(file, 1)
    statistic = load_data(file, 2)

    snimage = data.copy()
    snimage.data = data.data / statistic.data
    snimage.write("image00.fits")


def write_to_file(file_path, data):
    """
    Function to write data to a file

    Parameters:
    file_path (str): Path to the file
    data (float): Data to be written to the file
    """
    with open(file_path, "ab") as fout:
        fout.write(struct.pack("f", data))


def process_cube_data(cube_data, dy, dx, file_path):
    """
    Function to process cube data and write it to a file

    Parameters:
    cube_data (numpy array): The cube data to be processed
    dy (int): The y dimension of the cube data
    dx (int): The x dimension of the cube data
    file_path (str): Path to the file where the processed data will be written
    """
    with open(file_path, "ab") as fout:
        for y in range(dy):
            for x in range(dx):
                spec = cube_data[:, y, x]
                myfmt = "f" * len(spec.data)
                fout.write(struct.pack(myfmt, *spec.data))


create_masking_plot()

filename_second_argument = sys.argv[2]
processed_file = os.path.join(project_path_config.DATA_PATH_PROCESSED, filename_second_argument)

c0 = mpdaf.obj.Cube(processed_file)
dimension_z, dimension_y, dimension_x = c0.shape
start = c0.wave.get_crval()

write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), dimension_z)
write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), dimension_y)
write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), dimension_x)
write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), start)

process_cube_data(c0.data, dimension_y, dimension_x,
                  os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"))