import os
import sys
import project_path_config
import mpdaf


def create_masking_plot():
	filename = sys.argv[1]
	file = os.path.join(project_path_config.DATA_PATH_RAW, filename)

	data = mpdaf.obj.Cube(file, ext=1).sum(axis=0)
	data.info()
	print(f"data loaded\n")

	statistic = mpdaf.obj.Cube(file, ext=2).sum(axis=0)
	statistic.info()
	print(f"statistic loaded\n")

	snimage = data.copy()
	snimage.data = data.data / statistic.data
	snimage.info()
	snimage.write("image00.fits")
	print(f"file created")


create_masking_plot()
