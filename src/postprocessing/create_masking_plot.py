import os
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
    filename = sys.argv[1]
    file = os.path.join(project_path_config.DATA_PATH_RAW, filename)

    data = load_data(file, 1)
    statistic = load_data(file, 2)

    snimage = data.copy()
    snimage.data = data.data / statistic.data
    snimage.write("image00.fits")


create_masking_plot()
