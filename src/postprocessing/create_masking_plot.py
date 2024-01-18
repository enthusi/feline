import os
import sys
import mpdaf
import project_path_config


def load_data(file, ext):
    return mpdaf.obj.Cube(file, ext=ext).sum(axis=0)


def create_masking_plot():
    filename = sys.argv[1]
    file = os.path.join(project_path_config.DATA_PATH_RAW, filename)

    data = load_data(file, 1)
    statistic = load_data(file, 2)

    snimage = data.copy()
    snimage.data = data.data / statistic.data
    snimage.write("image00.fits")


create_masking_plot()
