import os
import struct
import sys
import mpdaf.obj
import project_path_config


def write_to_file(file_path, data):
    """
    Function to write data to a file

    Parameters:
    file_path (str): Path to the file
    data (float): Data to be written to the file
    """
    with open(file_path, "wb") as fout:
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
    with open(file_path, "wb") as fout:
        for y in range(dy):
            for x in range(dx):
                spec = cube_data[:, y, x]
                myfmt = "f" * len(spec.data)
                fout.write(struct.pack(myfmt, *spec.data))


filename = sys.argv[1]
file = os.path.join(project_path_config.DATA_PATH_PROCESSED, filename)

c0 = mpdaf.obj.Cube(file)
dz, dy, dx = c0.shape
start = c0.wave.get_crval()

write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), dz)
write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), dy)
write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), dx)
write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), start)

process_cube_data(c0.data, dy, dx, os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"))
