"""
This module, `masking_and_transpose`, is used for preprocessing and analyzing
astronomical data. It includes functions for data loading, masking, and cube
processing to prepare the data for further analysis.

Functions:
    - load_data: Loads and sums data from a file
      along a specified axis.
    - create_masking_plot: Generates a masking
      plot based on signal-to-noise ratio.
    - write_to_file: Writes numerical data to a file in binary
      format.
    - transpose_cube_data: transposes and writes spectra
      for each pixel for analysis.

Dependencies:
    - mpdaf: For handling astronomical data cubes.
    - numpy: For numerical calculations.
    - os, sys, struct: For file and system operations.
"""

import os
import struct
import sys
import numpy as np
import mpdaf

import src.project_path_config as project_path_config


def load_data(file: str, ext: int) -> mpdaf.obj.Cube:
    """
    Load data from a file.

    Args:
        file (str): The file path.
        ext (int): The extension number.

    Returns:
        mpdaf.obj.Cube: The sum of the cube data along the specified axis.
    """
    return mpdaf.obj.Cube(file, ext=ext).sum(axis=0)


def create_masking_plot() -> None:
    """
    Create a masking plot and write it to a file.

    The function takes no parameters.
    It uses the filename from the command line arguments,
    loads the data and statistic from the file,
    creates a signal-to-noise image, and writes it to 'image00.fits'.

    Returns:
        None
    """
    filename_first_argument = sys.argv[1]
    file = os.path.join(
        project_path_config.DATA_PATH_RAW,
        filename_first_argument)
    data = load_data(file, 1)
    statistic = load_data(file, 2)
    snimage = data.copy()
    snimage.data = data.data / statistic.data
    snimage.write(f"{project_path_config.DATA_PATH_PROCESSED}/image00.fits")


def write_to_file(file_path: str, data: float) -> None:
    """
    Function to write data to a file.

    Args:
        file_path (str): Path to the file
        data (float): Data to be written to the file
    Returns:
        None
    """
    with open(file_path, "ab") as fout:
        fout.write(struct.pack("f", data))


def transpose_cube_data(cube_data: np.ndarray, dy: int, dx: int,
                        file_path: str) -> None:
    """
    Function to transpose cube data and write it to a file.

    Args:
        cube_data (np.ndarray): The cube data to be processed
        dy (int): The y dimension of the cube data
        dx (int): The x dimension of the cube data
        file_path (str): Path to the file where the processed data will
                         be written
    Returns:
        None
    """
    with open(file_path, "ab") as fout:
        for y in range(dy):
            for x in range(dx):
                spec = cube_data[:, y, x]
                myfmt = "f" * len(spec.data)
                fout.write(struct.pack(myfmt, *spec.data))


if __name__ == "__main__":
    create_masking_plot()

    filename_second_argument = sys.argv[2]
    processed_file = os.path.join(
        project_path_config.DATA_PATH_PROCESSED,
        filename_second_argument)

    cube_from_processed_file = mpdaf.obj.Cube(processed_file, ext=1)
    cube_dimensions_z, cube_dimensions_y, cube_dimensions_x = (
        cube_from_processed_file.shape)
    wavelength_at_reference_pixel = cube_from_processed_file.wave.get_crval()

    write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED,
                               "raw_reordered_s2ncube.dat"),
                  cube_dimensions_z)
    write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED,
                               "raw_reordered_s2ncube.dat"),
                  cube_dimensions_y)
    write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED,
                               "raw_reordered_s2ncube.dat"),
                  cube_dimensions_x)
    write_to_file(os.path.join(project_path_config.DATA_PATH_PROCESSED,
                               "raw_reordered_s2ncube.dat"),
                  wavelength_at_reference_pixel)

    transpose_cube_data(cube_from_processed_file.data, cube_dimensions_y,
                        cube_dimensions_x,
                        os.path.join(project_path_config.DATA_PATH_PROCESSED,
                                   "raw_reordered_s2ncube.dat"))
