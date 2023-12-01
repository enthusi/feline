import sys
import mpdaf.obj

first = True

if first:
	name = sys.argv[1]
	data = mpdaf.obj.Cube(name, ext=1).sum(axis=0)
	data.info()
	print("data loaded")
	stat = mpdaf.obj.Cube(name, ext=2).sum(axis=0)
	stat.info()
	print("stat loaded")
	snimage = data.copy()
	snimage.data = data.data / stat.data
	snimage.info()
	snimage.write('image00.fits')
else:
	snimage = mpdaf.obj.Image('image00.fits')
