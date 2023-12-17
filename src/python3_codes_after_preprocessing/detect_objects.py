import os
import struct
import sys
import astropy
import imageio
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
import skimage

mpl.use('TkAgg')


def onclick(event):
	print(('#button=%d, x=%d, y=%d, xdata=%f, ydata=%f' % (event.button, event.x, event.y, event.xdata, event.ydata)))
	manuals = []  # was ist mit der manuals liste?
	plt.plot(event.xdata, event.ydata, 'yo')
	manuals.append((event.ydata, event.xdata))  # was ist mit der manuals liste?
	plt.show()


def world_to_pix(coord, rad):
	radarray = np.array([[rad[0], rad[1], 0]], np.float_)
	world = coord.wcs_world2pix(radarray, 0)
	x = world[0][0]
	y = world[0][1]
	return x, y


def pix_to_world(coord, pix):
	pixarray = np.array([[pix[0], pix[1], 0]], np.float_)
	world = coord.wcs_pix2world(pixarray, 0)
	ra = world[0][0]
	dec = world[0][1]
	return ra, dec


hdu = astropy.io.fits.open(sys.argv[1])
coord = astropy.wcs.WCS(hdu[0].header)

with open("raw_reordered_s2ncube.dat", 'rb') as f:
	header = f.read()[:16]

dz = struct.unpack('f', header[0:4])[0]
xd = struct.unpack('f', header[4:8])[0]
yd = struct.unpack('f', header[8:12])[0]

xd = int(xd)
yd = int(yd)
dz = int(dz)

print("#Cube dimensions (z,y,x): %d, %d, %d" % (dz, xd, yd))

atoms = [
	[6564.61],
	[4862.72],
	[4341.68],
	[4102.89],
	[3727.09, 3729.88],
	[4960.30, 5008.24],
	[6549.86, 6585.27],
	[6718.29, 6732.67],
	[3869.81, 3968.53],
	[1908.10, 1906.05],
	[1215.67]
]

atom_id = {
	6564.61: "Ha",
	4862.72: "Hb",
	4341.68: "Hg",
	4102.89: "Hd",
	3727.09: "OIIa",
	3729.88: "OIIb",
	4960.30: "[OIII]a",
	5008.24: "[OIII]b",
	6549.86: "[NII]a",
	6585.27: "[NII]b",
	6718.29: "[SII]a",
	6732.67: "[SII]b",
	3869.81: "NeIIIa",
	3968.53: "NeIIIb",
	1908.10: "CIIIa",
	1906.05: "CIIIb",
	1215.45: "Lya"}

mpl.rcParams["savefig.directory"] = "."


def gauss2d(xy, amp, x0, y0, a, b, c):
	x, y = xy
	inner = a * (x - x0) ** 2
	inner += 2 * b * (x - x0) ** 2 * (y - y0) ** 2
	inner += c * (y - y0) ** 2
	return amp * np.exp(-inner)


def twoD_Gaussian(xxx_todo_changeme, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
	(x, y) = xxx_todo_changeme
	xo = float(xo)
	yo = float(yo)
	a = (np.cos(theta) ** 2) / (2 * sigma_x ** 2) + (np.sin(theta) ** 2) / (2 * sigma_y ** 2)
	b = -(np.sin(2 * theta)) / (4 * sigma_x ** 2) + (np.sin(2 * theta)) / (4 * sigma_y ** 2)
	c = (np.sin(theta) ** 2) / (2 * sigma_x ** 2) + (np.cos(theta) ** 2) / (2 * sigma_y ** 2)
	g = offset + amplitude * np.exp(- (a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2)))
	return g.ravel()


def print_lines(toggle, z):
	for k in range(len(atoms)):
		# is k in the template?
		if toggle & 0x1 == 0:
			toggle = toggle // 2
			continue

		# ok, we consider this atom/transition
		toggle = toggle // 2
		atom = atoms[k]

		for emission in atom:
			pos = emission * (z + 1)
			name = atom_id[emission]
			print("%s (%.1f)," % (name, pos), end=' ')
	print()


isize = xd * yd
size = isize

data = np.fromfile('float32_array_omp4.raw', dtype='float32')
plane, redshift, template, used = np.split(data, 4)

plane.resize((xd, yd))
redshift.resize((xd, yd))
template.resize((xd, yd))
used.resize((xd, yd))

# Scale the floating-point values to the range [0, 1]
plane_scaled = (plane - np.min(plane)) / (np.max(plane) - np.min(plane))

# Convert the floating-point values to integers in the range [0, 255]
plane_uint8 = (plane_scaled * 255).astype(np.uint8)

# Save the image using imageio.imwrite
imageio.imsave('image.png', plane_uint8)

if os.path.isfile("imagemask.png"):
	mymask = imageio.v2.imread("imagemask.png")

	ny, nx = mymask.shape
	for iy in range(ny):
		for ix in range(nx):
			if mymask[iy, ix] == 0xff:
				plane[iy, ix] = 0

else:
	pass
# print 'CAUTION! no valid imagemask.png defined yet!'

data = scipy.ndimage.gaussian_filter(plane, sigma=1)

# data=plane
coordinates = skimage.feature.peak_local_max(data, min_distance=1, threshold_abs=50, exclude_border=10, num_peaks=300)

width = 0
bleft = width
btop = width
bright = xd - width
bbottom = yd - width

hdu_muse = astropy.io.fits.open("image00.fits", memmap=False)
data_muse = hdu_muse[1].data
nan_sel = np.isnan(data_muse)

xy = []
ay = []
ax = []
for val in coordinates:
	nx, ny = val
	ax.append(nx)
	ay.append(ny)
	xy.append((nx, ny))

print("#", len(xy))

test = plt.imshow(data, vmin=30, vmax=500, interpolation='nearest', cmap='jet')

plt.colorbar()

plt.autoscale(False)

run_id = 0
for hit in xy:  # y,x
	y, x = hit
	if x < 1: continue
	y = int(y)
	x = int(x)
	u_i = int(used[y, x])
	z_i = redshift[y, x]
	q_i = plane[y, x]
	t_i = int(template[y, x])

	print("%d %d %d %1.6f %d %d %d" % (run_id, int(y), int(x), z_i, q_i, u_i, t_i), end=' ')
	ra, dec = pix_to_world(coord, (x, y))
	print("\t%.6f %.6f" % (ra, dec), end=' ')
	print_lines(t_i, z_i)
	run_id += 1

plt.plot(ay, ax, 'rx', markersize=2)

plt.title("%d sources " % (len(xy)))
plt.show()
plt.savefig('result.png', bbox_inches='tight')
