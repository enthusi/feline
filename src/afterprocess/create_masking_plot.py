import sys
import mpdaf
from config import PROJECT_ROOT

ENABLE = True

if ENABLE:
	name = sys.argv[1]
	data = mpdaf.obj.Cube(name, ext=1).sum(axis=0)
	data.info()
	print("data loaded")
	statistic = mpdaf.obj.Cube(name, ext=2).sum(axis=0)
	statistic.info()
	print("statistic loaded")
	snimage = data.copy()
	snimage.data = data.data / statistic.data
	snimage.info()
	snimage.write('image00.fits')
else:
	snimage = mpdaf.obj.Image('image00.fits')
