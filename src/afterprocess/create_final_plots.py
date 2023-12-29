import mpdaf
import matplotlib.pyplot as plt
import math
import struct
import numpy as np
import scipy
import ref_index
import matplotlib as mpl
import sys
import astropy.cosmology
import astropy.io.fits
import astropy.wcs
import project_path_config
import os

# from mpdaf.obj import Cube
# from numpy import *
# from scipy.optimize import curve_fit
# from astropy import wcs
# from astropy.io import fits

mpl.use('TkAgg')

# from mpdaf.obj import CubeDisk
# this is the version as of 2017-11-06!
# python new_plot_and_fit_catalog_allz.py cube.fits s2n_v250.fits sorted_catalog.txt med_filt.fits HE0153-4520

cosmo = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

# qsotag=sys.argv[5].ljust(14)
galaxy_absorption = [
	2586.65,
	2600.17,
	2796.35,
	3890.1506,
	3934.777,
	3969.588,
	4102.89,
	4305.61,
	4341.68,
	4862.68,
	5176.7,
	5895.6
]

"""
| Transition  |  wl (vac) | comment                                 |
|-------------+-----------+-----------------------------------------|
| Fe II       |   2586.65 |                                         |
| Fe II       |   2600.17 |                                         |
| Mg II       |   2796.35 | no need to mark 2803, I think           |
| H$\zeta$    | 3890.1506 |                                         |
| K(Ca)       |  3934.777 |                                         |
| H(Ca)       |  3969.588 |                                         |
| H$\epsilon$ | 3971.1951 | blended with H(Ca), so maybe don't mark |
| H$\delta$   |   4102.89 |                                         |
| G           |   4305.61 |                                         |
| $H\gamma$   |   4341.68 |                                         |
| $H\beta$    |   4862.68 |                                         |
| MgI (blend) |    5176.7 |                                         |
| Na          |    5895.6 |  
"""

quiet = False


def scale_params(redshift):
	ps = cosmo.kpc_proper_per_arcmin(redshift).value / 60.0
	return ps


# that function gives you the plate scale at the redshift of the galaxy
# to convert arcsec to kpc
# then:

def get_impact(QSO_X, QSO_Y, px, py, z):
	print(QSO_X, QSO_Y, px, py, z)
	theta = math.sqrt((QSO_X - px) ** 2 + (QSO_Y - py) ** 2) * 0.2
	scale = scale_params(z)
	# print theta,scale, theta*scale, b
	return theta * scale


# dummies
qso_x = 150
qso_y = 150
foundb = -1
foundz = -1
foundm = -1
foundqop = -1
foundifrom = "not"

mpl.rcParams["savefig.directory"] = "."
cont_flag = False

full_data = open('info.txt', 'w')
full_data.write('id, px, py, z, b, quality, template, used, gaussx, gaussy,' ', gaussrot, smoothratio, oiiratio, gal_ra, gal_dec, qso_x, qso_y\n')

error_log = open('errors.log', 'w')


def pix_to_world(coord, pix):
	print(pix)

	pixarray = np.array([[pix[0], pix[1], 0]], np.float_)
	print(pixarray)
	world = coord.wcs_pix2world(pixarray, 0)
	ra = world[0][0]
	dec = world[0][1]
	return ra, dec


def world_to_pix(coord, rad):
	# print pix

	radarray = np.array([[rad[0], rad[1], 0]], np.float_)
	# print pixarray
	world = coord.wcs_world2pix(radarray, 0)
	x = world[0][0]
	y = world[0][1]
	return x, y


global px, py, z, run_id, forfit_t, forfit_w, used, gtemplate, npx, npy, zerr, ra, dec, verr, quality
forfit_t = 0
forfit_w = np.zeros(10)
global use_new_pos
use_new_pos = 0

prev_cat = False

# 12 just woudlnt fit in!
max_lines_shown = 12
columns = 12
rows = 4
ratios = open(project_path_config.DATA_PATH_LOOKUP + "stat_oii_ratio.txt", "w")

with open(project_path_config.DATA_PATH_PROCESSED + "raw_reordered_s2ncube.dat", 'rb') as f:
	header = f.read()[:16]

dz = struct.unpack('f', header[0:4])[0]
xd = struct.unpack('f', header[4:8])[0]
yd = struct.unpack('f', header[8:12])[0]
crval = struct.unpack('f', header[12:16])[0]
crmax = crval + dz * 1.25

dz = int(dz)
xd = int(xd)
yd = int(yd)

print("#Cube dimensions (z,y,x): %d, %d, %d" % (dz, xd, yd))

Lorrie = False

read_own = False
min_line_number = 2

only_good = False

statistics = open('statistics.txt', 'w')
peakratio = 0

try:
	pcatalog = {}
	previous = open('input_catalog.txt', 'r')
	for line in previous:
		run_id = int((line.split()[0]))
		mark = int((line.split()[4]))
		pcatalog[run_id] = mark
	prev_cat = True
	previous.close()
	print("read in former selections")
except:
	print("no former selections available!")

verified_sources = open('output_catalog.txt', 'w')
export_ds9 = True


def fill(orig, start_coords, fill_value):
	area = 0
	"""
	Flood fill algorithm
	
	Parameters
	----------
	data : (M, N) ndarray of uint8 type
		Image with flood to be filled. Modified inplace.
	start_coords : tuple
		Length-2 tuple of ints defining (row, col) start coordinates.
	fill_value : int
		Value the flooded area will take after the fill.
		
	Returns
	-------
	None, ``data`` is modified inplace.
	"""
	data = orig.copy()
	xsize, ysize = data.shape
	orig_value = data[start_coords[0], start_coords[1]]

	stack = {(start_coords[0], start_coords[1])}
	if fill_value == orig_value:
		print("Filling region with same value ")
		return 0

	while stack:
		x, y = stack.pop()

		if data[x, y] == orig_value:
			data[x, y] = fill_value
			area += 1
			if x > 0:
				stack.add((x - 1, y))
			if x < (xsize - 1):
				stack.add((x + 1, y))
			if y > 0:
				stack.add((x, y - 1))
			if y < (ysize - 1):
				stack.add((x, y + 1))
	return area


def get_num_lines(toggle):
	lines = 0
	for k in range(len(atoms)):
		if toggle & 0x1 == 0:
			toggle = toggle / 2
			continue
		toggle = toggle / 2
		atom = atoms[k]
		# atoms_found.append(k)
		for emission in atom:
			lines += 1
	return lines


def gauss_function(x, a, x0, sigma):
	return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def galaxy(w, *p):
	global forfit_t, atoms
	z = p[0]
	sigma = p[1]
	# print p
	toggle = forfit_t
	flux = np.zeros(len(w))
	i = 0
	for k in range(len(atoms)):
		if toggle & 0x1 == 0:
			toggle = toggle / 2
			continue
		toggle = toggle / 2
		atom = atoms[k]
		# atoms_found.append(k)
		for emission in atom:
			vacline = emission
			pos = vacline * (z + 1)
			newpos = ref_index.vac2air(pos / 10.0) * 10.0

			amplitude = p[2 + i]

			flux += gauss_function(w, amplitude, newpos, sigma)

			i += 1
	return flux


def fit_template(t, z, f, w, sigma_array, scipy):
	global forfit_t, forfit_w
	forfit_t = t
	forfit_w = w

	params = []

	params.append(z)
	params.append(1.0)
	param_bounds_low = []
	param_bounds_high = []
	# how many atoms?
	# count in t
	param_bounds_low.append(z - 0.002)
	param_bounds_high.append(z + 0.002)

	param_bounds_low.append(0.9)
	param_bounds_high.append(4.0)

	# but how many actual lines are that?
	lines = get_num_lines(t)

	for i in range(lines):
		amp = 20.0
		sig = 1.0
		params.append(amp)

		param_bounds_low.append(0)  # amp
		param_bounds_high.append(np.inf)

	param_bounds = (param_bounds_low, param_bounds_high)
	popt, pcov = scipy.optimize.curve_fit(galaxy, w, f, p0=params, bounds=param_bounds, max_nfev=1000)

	# print popt,pcov
	try:
		perr = np.sqrt(np.diag(pcov))[0]
	except:
		perr = 80

	new_z = popt[0]
	print("------------------------")
	print(z, new_z, (new_z - z) * 300000, t, perr)
	# print popt
	print(popt)
	return new_z, perr, popt


def b0pressed(event):
	correct_pos()
	print("ok", px, py, z)
	if foundz > 0:
		verified_sources.write(
			"%d \t%.1f \t%.1f \t%f  1 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%.1f\t%.1f\n" % (
				run_id, px, py, z, verr, used, gtemplate, found, ra, dec, verr, quality, z, foundz, foundm,
				get_impact(qso_x, qso_y, px, py, z), foundb))
		stat_good.write("%f %f %f\n" % (z, quality, peakratio))
	else:
		verified_sources.write("%d \t%.1f \t%.1f \t%f  1 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\t-\t-\t-\t%.1f\n" % (
			run_id, px, py, z, verr, used, gtemplate, found, ra, dec, verr, quality,
			get_impact(qso_x, qso_y, px, py, z)))
	verified_sources.flush()
	plt.close()


def b1pressed(event):
	print("bad:", px, py, z)
	verified_sources.write("%d \t%.1f \t%.1f \t%f  0 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\n" % (
		run_id, px, py, z, verr, used, gtemplate, found, ra, dec, verr, quality))
	verified_sources.flush()
	stat_bad.write("%f %f %f\n" % (z, quality, peakratio))
	plt.close()


def b2pressed(event):
	wave = ref_index.vac2air((1 + z) * 3727.09 / 10.0) * 10.0
	lyz = ref_index.air2vac(wave / 10.0) * 10.0 / 1215.67 - 1
	print("OIIa marked as Lyman alpha ", px, py, z, lyz)
	verified_sources.write("%d \t%.1f \t%.1f \t%f  2 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\n" % (
		run_id, px, py, lyz, verr, used, gtemplate, found, ra, dec, verr, quality))
	verified_sources.flush()
	plt.close()


def b2pressedb(event):
	# correct_pos()

	# compute a  bad new z assuming Lyman alpha
	# from z back to wave
	wave = ref_index.vac2air((1 + z) * 3729.88 / 10.0) * 10.0
	lyz = ref_index.air2vac(wave / 10.0) * 10.0 / 1215.67 - 1
	print("OIIb marked as Lyman alpha ", px, py, z, lyz)
	verified_sources.write("%d \t%.1f \t%.1f \t%f  2 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\n" % (
		run_id, px, py, lyz, verr, used, gtemplate, found, ra, dec, verr, quality))
	verified_sources.flush()
	plt.close()


def b3pressed(event):
	if pcatalog[run_id] > 0:  # dont correct positions for 'NO'
		correct_pos()
	print("use old", px, py, z)
	verified_sources.write("%d \t%.1f \t%.1f \t%f %d \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\n" % (
		run_id, px, py, z, pcatalog[run_id], verr, used, gtemplate, found, ra, dec, verr, quality))
	verified_sources.flush()
	plt.close()


def correct_pos():
	if read_own:
		print("wont fix own list")
		return

	global px, py, npx, npy
	use_new_pos = check.lines[0][0].get_visible()
	print("===>", use_new_pos)
	if use_new_pos > 0:
		px = npx
		py = npy
		ra, dec = pix_to_world(coord, (px, py))


def checkpressed(event):
	# global use_new_pos
	# if event == '2 Hz':
	#  use_new_pos^=1
	#  print use_new_pos
	pass


def correctlimit(ax, x, y):
	# ax: axes object handle
	#  x: data for entire x-axes
	#  y: data for entire y-axes
	# assumption: you have already set the x-limit as desired

	lims = ax.get_xlim()

	i = np.where((x > lims[0]) & (x < lims[1]))[0]
	range = y[i].max() - y[i].min()
	ax.set_ylim(-0.4 * range, 1.4 * range)


if len(sys.argv) < 2:
	print("SYNTAX: %s cube.fits catalog.cat [ds9.reg]" % sys.argv[0])
	sys.exit(0)

if len(sys.argv) == 6:
	regfile = open(sys.argv[5], 'w')
	size = 4.0
	export_ds9 = True
	regfile.write("""# Region file format: DS9 version 4.1
    # Filename: qso1422.fits[DATA]
    global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1
    image
    """)

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
	[1908.73, 1906.68],
	[1215.67]
]

atomsize = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1]

atom_id = {
	6564.61: r"H$\alpha$",
	4862.72: r"H$\beta$",
	4341.68: r"H$\gamma$",
	4102.89: r"H$\delta$",
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
	1908.73: "CIIIa",
	1906.68: "CIIIb",
	1215.67: "Lya"
}

data = np.fromfile(project_path_config.DATA_PATH_ROOT + "float32_array_omp4.raw", dtype="float32")
plane, redshift, template, imused = np.split(data, 4)

plane.resize((xd, yd))
redshift.resize((xd, yd))
template.resize((xd, yd))
imused.resize((xd, yd))

# for data cube
cube = mpdaf.obj.Cube(project_path_config.DATA_PATH_PROCESSED + sys.argv[4], ext=0)
cubestat = mpdaf.obj.Cube(project_path_config.DATA_PATH_PROCESSED + sys.argv[4], ext=1)
cube.info()

original_cube = mpdaf.obj.Cube(project_path_config.DATA_PATH_RAW + sys.argv[1], ext=1)
# NEW2017
# if WHITEimage is an extention
# whiteimage=Image(sys.argv[1],ext=4).data
# fullwhiteimage=Image(sys.argv[1],ext=4)
# use this for stand alone white image
whiteimage = mpdaf.obj.Cube(project_path_config.DATA_PATH_RAW + sys.argv[1], ext=1).sum(axis=0).data

# whiteimage=Image(sys.argv[1]).data
fullwhiteimage = mpdaf.obj.Cube(project_path_config.DATA_PATH_RAW + sys.argv[1], ext=1).sum(axis=0)

s2ncube = mpdaf.obj.Cube(project_path_config.DATA_PATH_PROCESSED + sys.argv[2], ext=0)
hdu = astropy.io.fits.open(project_path_config.DATA_PATH_PROCESSED + sys.argv[2])
coord = astropy.wcs.WCS(hdu[0].header)

dz, dy, dx = cube.shape

catalog = open(sys.argv[3])
colors = mpl.cm.get_cmap("winter")
colors._init()

i = 0
vp = 4
objid = 0
oldid = 0
objects = []
# spatial aperture (pix)
ds = 3

# interval blue and redward of detection (pix)
w = 15
try:
	qso_id = sys.argv[5]

except:
	print("no QSO given")
	sys.exit(1)

qso_positions_file = open(project_path_config.DATA_PATH_LOOKUP + "qso_centers_deg.txt", "r")
qso_found = False
print("looking for %s" % qso_id)
for qso in qso_positions_file:
	if qso_id in qso:
		qso_found = True
		qso_ra = float(qso.split()[1])
		qso_dec = float(qso.split()[2])

if not qso_found:
	print("no QSO found")
	sys.exit(1)

print(qso_ra)
print(qso_dec)
qso_x, qso_y = world_to_pix(coord, (qso_ra, qso_dec))
print(qso_x, qso_y)

for line in catalog:
	if line[0] == '#': continue
	# 496 294 25 1.017953 71 2 16 OIIa (7521.1), OIIb (7526.7),

	use_new_pos = 0

	# reset the values from Lorrie's catalog in case none is found!
	found = False
	foundb = -1
	foundz = -1
	foundm = -1
	foundqop = -1
	foundifrom = "not"
	# if reading in my own catalog
	if True:
		run_id = int((line.split()[0]))
		py = float((line.split()[1]))
		px = float((line.split()[2]))
		z = float(line.split()[3])

		quality = float((line.split()[4]))
		used = int((line.split()[5]))
		gtemplate = int((line.split()[6]))
		ra, dec = pix_to_world(coord, (px, py))
		# SPECIFIC for THIS cube
		border_distance = min(min(px, py), min(dx - px, dy - py))
		if border_distance < 15: continue

	if prev_cat:
		mark = pcatalog[run_id]

	# use THIS line to look for a specific object by ID
	# if quality > 300: continue
	# if min(px,py)<8:
	#  print "TOO CLOSE TO BORDER"
	#  continue

	# if run_id != 81: continue
	# start at a certain ID
	# if run_id==448:
	#    cont_flag = True
	#    continue
	# if not cont_flag: continue
	# if run_id in [18]:
	#    continue
	print()
	print("running id", run_id)
	if only_good:
		if mark == 0: continue
	# if run_id==289: continue
	# if run_id==216: continue
	# if run_id==450: continue

	if used < min_line_number: continue
	# if quality < 200:continue
	toggle = gtemplate
	positions = []

	atoms_found = []
	lines_found = []

	raw_flux = cube[:, int(py) - ds:int(py) + ds, int(px) - ds:int(px) + ds].mean(axis=(1, 2))
	raw_data = raw_flux.data
	# raw_wave=np.arange(raw_flux.wave.crval, raw_flux.wave.crval+raw_flux.wave.cdelt*raw_flux.wave.shape,raw_flux.wave.cdelt)
	raw_wave = np.arange(raw_flux.wave.get_crval(),
						 raw_flux.wave.get_crval() + raw_flux.wave.get_step() * raw_flux.wave.shape,
						 raw_flux.wave.get_step())
	print("fitting now")
	raw_sigma = cubestat[:, int(py) - ds:int(py) + ds, int(px) - ds:int(px) + ds].mean(axis=(1, 2))
	# raw_sigma.data=sqrt(raw_sigma.data)#/math.sqrt((4*ds*ds)))
	# raw_flux.plot()
	# raw_sigma.plot()
	# plt.show()
	# raw_input("tets")

	# print raw_data/raw_sigma
	# sys.exit(1)

	valid_model = True
	try:

		newz, zerr, gal_model = fit_template(gtemplate, z, raw_data, raw_wave, raw_sigma.data)
	# print 1/0
	except:
		print("** fit did not converge!")
		error_log.write('no valid model %d %.3f %d\n' % (run_id, z, gtemplate))
		# sys.exit(1)
		newz = z
		zerr = 100
		gal_model = 0
		valid_model = False

	for k in range(len(atoms)):
		# is k in the template?
		if toggle & 0x1 == 0:
			toggle = toggle // 2
			continue

		# ok, we consider this atom/transition
		toggle = toggle // 2
		atom = atoms[k]
		atoms_found.append(k)
		for emission in atom:
			lines_found.append(emission)
			pos = emission * (z + 1)
			name = atom_id[emission]
			# print "%s (%.1f)," % (name,pos),
			positions.append(pos)
	# print
	print(positions)

	if export_ds9 and objid > oldid:
		regfile.write("circle (%d,%d,%f) # text = {%s}\n" % (px, py, size, objid))
	oldid = objid

	wavemin = 4780
	wavemax = 9300

	print("new plot", i)
	j = 0
	count = len(positions)
	c = 299792.458
	dv = (z - newz) * c
	verr = zerr * c / (z + 1)
	zguess = z
	z = newz

	plt.figure(figsize=(16, 9))

	ax1 = plt.subplot2grid((rows, columns), (0, 0), colspan=9)
	plt.title("%s id=%04d, x=%.1f, y=%.1f, ra=%.6f dec=%.6f z=%.6f" % (qso_id, run_id, px, py, ra, dec, z))
	# print "***",z,new_z,perr
	# print "delta: %.1f error on new: %.1f" % ((z-new_z)*c,perr*c)

	# plt.title("id=%d, x=%.1f, y=%.1f, z=%.6f (verr=%.1f), used=%d, quality=%d" % (run_id,px,py,z,verr,used,quality))

	ax2 = plt.subplot2grid((rows, columns), (1, 0), colspan=9)
	plt.title("(verr=%.1f), %d used lines, match strength=%d, b=%.1f" % (
		verr, used, quality, get_impact(qso_x, qso_y, px, py, z)))

	ax2.tick_params(
		axis='both',  # changes apply to the x-axis
		which='both',  # both major and minor ticks are affected
		bottom='on',  # ticks along the bottom edge are off
		top='off',  # ticks along the top edge are off
		labelbottom='off',
		right='off',
		left='on',
		labelleft='on')  # labels along the bottom edge are off

	# plot regions of absorption first
	for absline in galaxy_absorption:
		abs_wav = ref_index.vac2air(absline * (z + 1) / 10.0) * 10.0
		if abs_wav > crval and abs_wav < crmax:
			ax1.axvline(x=abs_wav, color='aquamarine', linestyle='-', linewidth=4.0)

	for absline in galaxy_absorption:
		abs_wav = ref_index.vac2air(absline * (z + 1) / 10.0) * 10.0
		if abs_wav > crval and abs_wav < crmax:
			ax2.axvline(x=abs_wav, color='aquamarine', linestyle='-', linewidth=4.0)

	# plot actual flux spectrum
	spec = cube[:, int(py) - ds:int(py) + ds, int(px) - ds:int(px) + ds].mean(axis=(1, 2))
	original_spec = original_cube[:, int(py) - ds:int(py) + ds, int(px) - ds:int(px) + ds].mean(axis=(1, 2))
	# ax1=fig.add_subplot(gs[0])
	data1 = spec.data
	original_data1 = original_spec.data
	waven = np.arange(spec.wave.get_crval(), spec.wave.get_crval() + spec.wave.get_step() * spec.wave.shape,
					  spec.wave.get_step())
	ax1.step(waven, original_data1, where='mid', color='darkgrey')

	# plot possible positions for emission

	for a in atoms:
		for b in a:
			p = ref_index.vac2air(b * (z + 1) / 10.0) * 10.0
			if p > crval and p < crmax:
				ax1.axvline(x=p, color='k', linestyle='--')

	# plot detected emission  above that
	for g in range(len(positions)):
		wave = lines_found[g]

		thision = atom_id[wave]
		print(thision)
		wobs = ref_index.vac2air(wave * (z + 1) / 10.0) * 10.0
		print(wobs, positions[g])
		# wobs=ref_index.vac2air(wobs/10.0)*10.0
		ax1.axvline(x=wobs, color='r', linestyle='--')

	# plot possible positions for emission
	for a in atoms:
		for b in a:
			# p=b*(z+1)
			p = ref_index.vac2air(b * (z + 1) / 10.0) * 10.0
			if p > crval and p < crmax:
				ax2.axvline(x=p, color='k', linestyle='--')

	# plot detected emission  above that
	for g in range(len(positions)):
		wave = lines_found[g]
		thision = atom_id[wave]
		print(thision)
		wobs = ref_index.vac2air(wave * (z + 1) / 10.0) * 10.0
		print(wobs, positions[g])
		# gen narrow bands
		f = 4  # narrow band width
		d = 10
		p = int((wobs - crval) / 1.25)
		# all_band = cube[p-f:p+f, py-d:py+d ,px-d:px+d]
		all_band = cube[p - f:p + f, :, :]
		if g == 0:
			all_ima = all_band.sum(axis=0)
		else:
			all_ima += all_band.sum(axis=0)

		# wobs=ref_index.vac2air(wobs/10.0)*10.0
		ax2.axvline(x=wobs, color='r', linestyle='--')

	# var1=spec.var
	# print spec.var
	# sys.exit(1)

	waven_high = np.arange(spec.wave.get_crval(), spec.wave.get_crval() + spec.wave.get_step() * spec.wave.shape,
						   spec.wave.get_step() / 10.0)
	print("xxx")

	# waven=np.arange(spec.wave.crval, spec.wave.crval+spec.wave.cdelt*spec.wave.shape,spec.wave.cdelt)

	# ax1.step(waven,data1,where='mid')
	# ax1.set_xticks(np.arange(4699.59.0,9300.0,200))
	# waven_high=np.arange(spec.wave.crval, spec.wave.crval#+spec.wave.cdelt*spec.wave.shape,spec.wave.cdelt/10.0)
	# print "xxx"
	print(gal_model)
	print(forfit_t)
	print("xxx")

	ax1.step(waven, data1, where='mid', color='blue')
	try:
		ax1.plot(waven_high, galaxy(waven_high, *gal_model), 'k-')
	except:
		print("no model to plot available")
	ax1.set_xticks(np.arange(crval, crmax, 200))
	ax1.set_xlim(crval, crmax)
	# find bottom:
	lowest = min(data1)
	if lowest >= 0: bottom = -10
	if lowest < 0: bottom = lowest * 1.2

	ax1.set_ylim(bottom, max(data1) * 1.2)

	s2nspec = s2ncube[:, int(py) - ds:int(py) + ds, int(px) - ds:int(px) + ds].mean(axis=(1, 2))
	# test=fig.add_subplot(gs[1])
	data2 = s2nspec.data
	# waven=np.arange(spec.wave.crval, spec.wave.crval+spec.wave.cdelt*spec.wave.shape,spec.wave.cdelt)
	ax2.step(waven, data2, where='mid')
	ax2.set_xticks(np.arange(crval, crmax, 200))
	ax2.set_xlim(crval, crmax)

	# find bottom:
	lowest = min(data2)
	if lowest >= 0: bottom = -10
	if lowest < 0: bottom = lowest * 1.2

	ax2.set_ylim(bottom, max(data2) * 1.2)

	lines_found.sort()

	# plot all found lines AND always Ha,Hb
	hain = False
	hbin = False
	first = True
	height2b = 1
	height2a = 1

	oiifound = False
	for h in range(min(len(positions), max_lines_shown)):
		ax3 = plt.subplot2grid((rows, columns), (2, h))

		plt.title("%s" % (atom_id[lines_found[h]]), fontsize=10)
		wave = lines_found[h]
		# remember Oiii ratio

		#  oxy1=
		thision = atom_id[wave]
		print(thision)
		wobs = ref_index.vac2air(wave * (z + 1) / 10.0) * 10.0
		# wobs2=ref_index.vac2air(wave*(zguess+1)/10.0)*10.0
		print(wobs, positions[h])
		# wobs=ref_index.vac2air(wobs/10.0)*10.0
		ax3.axvline(x=wobs, color='r', linestyle='--')
		# ax3.axvline(x=wobs2,color='b', linestyle='--')

		# ax3.step(waven,data2,where='mid')
		try:
			ax3.plot(waven_high, galaxy(waven_high, *gal_model), 'r-')
		except:
			print("no model to plot available")
		ax3.plot(waven, data1, linestyle='-', drawstyle='steps-mid')
		ax3.plot(waven, data2, linestyle='-', drawstyle='steps-mid')
		fakewav = np.arange(wobs - 5, wobs + 5, 0.1)

		if valid_model:
			if atom_id[lines_found[h]] == "[OIII]a":
				heighta = sum(galaxy(fakewav, *gal_model))
			if atom_id[lines_found[h]] == "[OIII]b":
				heightb = sum(galaxy(fakewav, *gal_model))
				print("OOO", heightb / heighta, heighta, heightb)
				if heightb / heighta < 10: plt.xlabel("ratio:%.1f" % (heightb / heighta))
		else:
			heighta = 1
			heightb = 1

		if valid_model:
			if atom_id[lines_found[h]] == "OIIa":
				oiifound = True
				height2a = (galaxy([wobs], *gal_model))
			if atom_id[lines_found[h]] == "OIIb":
				height2b = (galaxy([wobs], *gal_model))

				ratios.write("%f\n" % (height2b / height2a))
				if height2b / height2a < 10: plt.xlabel("ratio:%.1f" % (height2b / height2a))
				oiifound = True
		else:
			height2a = 1
			height2b = 1

		if atom_id[lines_found[h]] == r"H$\alpha$": hain = True
		if atom_id[lines_found[h]] == r"H$\beta$": hbin = True
		dl = 15.0
		lim_low = max(crval, wobs - dl)
		lim_high = min(wobs + dl, crmax)
		while (lim_high - lim_low) < (2 * dl):
			lim_high += dl / 3.0
		ax3.set_xlim(lim_low, lim_high)
		ax3.set_xticks([wobs])

		# only plot axis/label for the first window
		if first:
			ax3.tick_params(
				axis='both',  # changes apply to the x-axis
				which='both',  # both major and minor ticks are affected
				bottom='on',  # ticks along the bottom edge are off
				top='off',  # ticks along the top edge are off
				labelbottom='on',
				right='off',
				left='on',
				labelleft='on')  # labels along the bottom edge are off

		else:
			ax3.tick_params(
				axis='both',  # changes apply to the x-axis
				which='both',  # both major and minor ticks are affected
				bottom='on',  # ticks along the bottom edge are off
				top='off',  # ticks along the top edge are off
				labelbottom='on',
				right='off',
				left='off',
				labelleft='off')  # labels along the bottom edge are off

		try:
			correctlimit(ax3, waven, data1)
		except:
			pass
		first = False

	# =================================================================
	# add Ha
	if False:  # hain==False and z<0.416:
		ax3 = plt.subplot2grid((rows, columns), (2, h + 1))
		plt.title("(Ha)", fontsize=10)
		wave = 6564.61
		# remember Oiii ratio

		#  oxy1=
		thision = atom_id[wave]
		print(thision)
		wobs = (wave * (z + 1))
		wobs2 = (wave * (zguess + 1))
		print(wobs, positions[h])
		# wobs=ref_index.vac2air(wobs/10.0)*10.0
		ax3.axvline(x=wobs, color='r', linestyle='--')
		# ax3.axvline(x=wobs2,color='b', linestyle='--')

		ax3.plot(waven, data1, linestyle='-', drawstyle='steps-mid')
		ax3.plot(waven, data2, linestyle='-', drawstyle='steps-mid')
		# ax3.step(waven,data2,where='mid')
		dl = 15.0
		ax3.set_xlim(wobs - dl, wobs + dl)
		ax3.set_xticks([wobs])
		correctlimit(ax3, waven, data1)
	# =================================================================
	plt.tight_layout()

	mark = 0
	if prev_cat:
		mark = pcatalog[run_id]

	# 2020 fix below ================================================
	galfit_x = 0
	galfit_y = 0
	galfit_rot = 0

	# NEW 2020! plot o2 zoom in with model
	if oiifound:
		ax4 = plt.subplot2grid((rows, columns), (3, 0), colspan=4)
		wobs = ref_index.vac2air(3728.0 * (z + 1) / 10.0) * 10.0
		# print wobs,positions[h]

		wobs1 = ref_index.vac2air(3727.09 * (z + 1) / 10.0) * 10.0
		wobs2 = ref_index.vac2air(3729.88 * (z + 1) / 10.0) * 10.0
		ax4.axvline(x=wobs1, color='k', linestyle='--')
		ax4.axvline(x=wobs2, color='k', linestyle='--')
		ax4.axhline(y=0, color='lightgrey')
		# ax3.axvline(x=wobs2,color='b', linestyle='--')

		# ax4.plot(waven,data2,linestyle='-', drawstyle='steps-mid')

		# ax3.step(waven,data2,where='mid')
		try:
			ax4.plot(waven_high, galaxy(waven_high, *gal_model), 'r-')

		except:
			print("no model to plot available")
		ax4.fill_between(waven, data1 - raw_sigma, data1 + raw_sigma, alpha=0.3, facecolor='#888888')
		ax4.plot(waven, data1, linestyle='-', drawstyle='steps-mid')
		# ax4.plot(waven,var1,linestyle='-', drawstyle='steps-mid')

		dl = 25.0
		lim_low = max(crval, wobs - dl)
		lim_high = min(wobs + dl, crmax)
		while (lim_high - lim_low) < (2 * dl):
			lim_high += dl / 3.0
		ax4.set_xlim(lim_low, lim_high)
		try:
			correctlimit(ax4, waven, data1)
		except:
			pass

		ax4.tick_params(
			axis='both',  # changes apply to the x-axis
			which='both',  # both major and minor ticks are affected
			bottom='on',  # ticks along the bottom edge are off
			top='off',  # ticks along the top edge are off
			labelbottom='off',
			right='off',
			left='off',
			labelleft='off')  # labels along the bottom edge are off

		try:
			correctlimit(ax4, waven, data1)
		except:
			pass

	# ================================
	if False:
		ax4 = plt.subplot2grid((rows, columns), (3, 0), colspan=2, rowspan=2)
		fluxes = []
		wcs1 = fullwhiteimage.wcs
		# for dc in range(3,20):
		dc = 10

		try:
			center_piece = Image(data=all_ima.data, wcs=wcs1)[int(py) - dc:int(py) + dc, int(px) - dc:int(px) + dc]
			seg = center_piece.segment(minsize=10)
			print(seg)
			source = seg[0]
			# source.plot(colorbar='v')
			# plt.show()
			# source.peak()
			# fwhmfit=source.fwhm()
			radius, ee = source.eer_curve(cont=source.background()[0])
			ax4.plot(radius, ee)
			plt.xlabel('radius')
			plt.ylabel('ERR')
			plt.tick_params(axis='both', left='off', top='off', right='off', bottom='on', labelleft='off',
							labeltop='off', labelright='off', labelbottom='on')
			# gfit = source.gauss_fit(maxiter=150, plot=True)
			# plt.tight_layout()
			del center_piece

			ax5 = plt.subplot2grid((rows, columns), (3, 2))
			gfit = source.gauss_fit(maxiter=150)
			# print "******************************"
			# try:
			#    print gfit.rot,gfit.fwhm ,"1"
			# except:
			#    print gfit.get_rot(),"2"
			gfitim = gauss_image(wcs=source.wcs, gauss=gfit)
			gresiduals = source - gfitim
			gfitim.plot()
			fwhmfit = gfit.fwhm

			plt.xlabel('x:%.2f, rot:%.1f' % (fwhmfit[1], gfit.rot))
			plt.ylabel('y:%.2f' % fwhmfit[0])
			galfit_x = fwhmfit[1]
			galfit_y = fwhmfit[0]
			galfit_rot = gfit.rot
			plt.title("2d gauss fit")
			plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off',
							labeltop='off', labelright='off', labelbottom='off')
			plt.gca().invert_yaxis()

			ax6 = plt.subplot2grid((rows, columns), (3, 3))
			# gfit = source.gauss_fit(maxiter=150)
			# gfitim = gauss_image(wcs=source.wcs, gauss=gfit)
			# gresiduals = source-gfitim
			gresiduals.plot()
			plt.title("gaussfit res")
			plt.xlabel('')
			plt.ylabel('')
			plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off',
							labeltop='off', labelright='off', labelbottom='off')
			plt.gca().invert_yaxis()
		# plt.tight_layout()
		except:
			error_log.write('could not segment source for %d\n' % run_id)
			print('could not segment source for %d' % run_id)
	# 2020 fix above ================================================

	bigpic = plt.subplot2grid((rows, 10), (0, 8), colspan=4, rowspan=2)
	bigpic.imshow(plane, vmax=1000, interpolation='none', cmap='jet')
	bigpic.plot(px, py, 'r*', ms=15)

	aw = 20
	aw = int(min(aw, px, py))

	wcs1 = fullwhiteimage.wcs
	narrowsa = mpdaf.obj.Image(data=all_ima.data, wcs=wcs1)[int(py) - aw // 2:int(py) + aw // 2,
			   int(px) - aw // 2:int(px) + aw // 2]
	spic = plt.subplot2grid((rows, columns), (3, 5))
	smoothnarrows = narrowsa.fftconvolve_gauss(center=None, flux=1.0, fwhm=(0.7, 0.7), peak=False, rot=0.0, factor=1,
											   unit_fwhm=None, inplace=False)
	# gaussian_filter(sigma=1,inplace=False)

	# narrowsa.write('testa.fits')
	# narrowsb.write('testb.fits')
	maxa = np.max(narrowsa.data)
	maxb = np.max(smoothnarrows.data)
	# ratio=narrowsa/smoothnarrows
	# ratio=narrowsa.correlate2d(smoothnarrows.data)
	# ratio=Image(data=ratio.data, wcs=wcs1)[int(py)-aw/2:int(py)+aw/2,int(px)-aw/2:int(px)+aw/2]
	# ratio.write('testc.fits')
	print("-=-=-=-=-=-=")
	# narrowsa.info()
	# narrowsb.info()
	# ratio.info()
	print("-=-=-=-=-=-=")

	# the following block is for the 2nd image but we need the peak first

	wcs1 = fullwhiteimage.wcs
	narrows = mpdaf.obj.Image(data=all_ima.data, wcs=wcs1)[int(py) - aw // 2:int(py) + aw // 2,
			  int(px) - aw // 2:int(px) + aw // 2]

	peakratio = float(maxa / maxb)

	# center_value=narrows[aw/2,aw/2]

	center_area = narrows[aw // 2 - 3:aw // 2 + 3, aw // 2 - 3:aw // 2 + 3]
	center_mean = np.mean(center_area.data)
	center_std = np.std(center_area.data)
	# center_peak=whitezoom[aw,aw]
	center_value = center_mean + center_std * 2.0

	plt.imshow(smoothnarrows.data, interpolation='none', cmap='jet', vmax=center_value)
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	plt.title("max %.1f" % (peakratio))

	spic = plt.subplot2grid((rows, columns), (3, 6))

	plt.imshow(narrows.data, interpolation='none', cmap='jet', vmax=center_value)
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	# plt.tight_layout()
	# for some mysterious reason this is required
	# plt.gca().invert_yaxis()
	# try:
	#  plt.colorbar()
	# except:
	#  pass
	# compute ratio of original divided by smoothed:

	plt.title("collapsed")

	# whiteimage

	whitezoom = whiteimage[int(py) - aw:int(py) + aw, int(px) - aw:int(px) + aw]
	spic = plt.subplot2grid((rows, columns), (3, 7))

	center_area = whitezoom[aw - 4:aw + 4, aw - 4:aw + 4]
	center_mean = np.mean(center_area.data)
	center_std = np.std(center_area.data)
	center_peak = whitezoom[aw, aw]
	center_value = center_mean + center_std * 4.0

	plt.imshow(whitezoom, interpolation='none', cmap='jet', vmax=center_value)
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	# plt.tight_layout()
	# try:
	#  plt.colorbar()
	# except:
	#  pass

	plt.title("white")
	plt.xlim(aw - 10, aw + 10)
	plt.ylim(aw - 10, aw + 10)
	# plt.ylim(py-10,py+10)
	# plt.axis('off')
	# plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off', labelright='off', labelbottom='off')
	# spic.plot(px,py,'r*',ms=10)
	# spic.plot(aw+px-int(px),aw+py-int(py),'r*',ms=10)
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	plt.gca().invert_yaxis()

	wcs1 = fullwhiteimage.wcs
	full_plane = mpdaf.obj.Image(data=plane, wcs=wcs1)[int(py) - aw:int(py) + aw, int(px) - aw:int(px) + aw]

	spic = plt.subplot2grid((rows, columns), (3, 8))
	# planbezoom=plane[py-aw:py+aw,px-aw:px+aw]
	center_value = full_plane[aw, aw]

	plt.imshow(full_plane.data, interpolation='none', cmap='jet', vmax=1.0 * center_value)
	# plt.tight_layout()

	# try:
	#  plt.colorbar()
	# except:
	#  pass
	# plt.gca().invert_yaxis()
	plt.title("quality")
	# plt.xlim(aw-10,aw+10)
	# plt.ylim(aw-10,aw+10)
	plt.xlim(aw - 10, aw + 10)
	plt.ylim(aw - 10, aw + 10)

	plt.axis('off')
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	# spic.plot(aw,aw,'r*',ms=10)
	# spic.plot(aw+px-int(px),aw+py-int(py),'r*',ms=10)
	plt.gca().invert_yaxis()
	# http://mpdaf.readthedocs.io/en/latest/api/mpdaf.obj.Gauss2D.html#mpdaf.obj.Gauss2D
	# bestgauss=full_plane.gauss_fit(circular=True,pix=True,pos_min=[aw-4,aw-4],pos_max=[aw+4,aw+10],plot=True)
	bestgauss = full_plane.gauss_fit(circular=True, pos_min=[aw - 4, aw - 4], pos_max=[aw + 4, aw + 10],
									 unit_center=None, unit_fwhm=None)
	a, b = bestgauss.center
	# spic.plot(aw+a-int(a),aw+b-int(b),'w.',ms=10)
	fwhm = bestgauss.fwhm

	# redshift
	spic = plt.subplot2grid((rows, columns), (3, 9))
	# plt.imshow(redshift, interpolation='none')
	# aw=10
	# aw=min(aw,px,py)
	testarea = redshift[int(py) - aw:int(py) + aw, int(px) - aw:int(px) + aw]

	plt.imshow(testarea, interpolation='none', cmap='jet')
	area = fill(testarea, (aw, aw), 2)
	plt.xlabel("area:%d" % area)
	# plt.xlim(px-10,px+10)
	# plt.ylim(py-10,py+10)
	plt.axis('on')
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	# spic.plot(px,py,'r*',ms=10)
	# spic.plot(aw,aw,'r*',ms=10)
	# try:
	#  plt.colorbar()
	# except:
	#  pass

	plt.title("redshift")

	# imused
	spic = plt.subplot2grid((rows, columns), (3, 10))
	plt.imshow(imused, interpolation='none', cmap='jet')

	plt.xlim(px - 10, px + 10)
	plt.ylim(py - 10, py + 10)
	plt.axis('off')
	plt.tick_params(axis='both', left='off', top='off', right='off', bottom='off', labelleft='off', labeltop='off',
					labelright='off', labelbottom='off')
	# spic.plot(px,py,'r*',ms=10)
	# try:
	#  plt.colorbar()
	# except:
	#  pass
	plt.gca().invert_yaxis()
	plt.title("no. lines")

	# go FULLSCREEN for the window(!)
	# this might only work in Linux!
	# http://stackoverflow.com/questions/12439588/how-to-maximize-a-plt-show-window-using-python
	# mng = plt.get_current_fig_manager()
	# mng.resize(*mng.window.maxsize())
	mng = plt.get_current_fig_manager()
	# mng.full_screen_toggle()

	print(px, py)
	npx = px - aw + b
	npy = py - aw + a
	print(npx, npy)

	# only write plots for verified ones!

	plt.savefig('%04d_%04d_%d_f%d_fig%04d.png' % (quality, run_id, mark, found, i))  # show()
	# plt.savefig('pdf%04d_%04d_%d_f%d_fig%04d.pdf'% (quality,run_id,mark,found,i))#show()
	i += 1
	if not quiet:
		plt.show()

	statistics.write(
		"%d %d %d %d %.1f %d %d %.1f\n" % (run_id, mark, area, border_distance, verr, used, use_new_pos, peakratio))
	statistics.flush()

	# compute RA,DEC from px,py
	# full_data.write('#id, px, py, z, b, quality, template, used, gaussx, gaussy, gaussrot, smoothratio, oiiratio, qso_ra, qso_dec, qso_x, qso_y')

	full_data.write('%04d, ' % run_id)
	full_data.write('%.2f, ' % px)
	full_data.write('%.2f, ' % py)
	full_data.write('%.6f, ' % z)
	full_data.write('%.1f, ' % get_impact(qso_x, qso_y, px, py, z))
	full_data.write('%d, ' % quality)
	full_data.write('%d, ' % gtemplate)
	full_data.write('%d, ' % used)

	full_data.write('%.2f, ' % galfit_x)
	full_data.write('%.2f, ' % galfit_y)
	full_data.write('%.2f, ' % galfit_rot)
	full_data.write('%.2f, ' % peakratio)
	full_data.write('%.2f, ' % (height2b / height2a))
	full_data.write('%.4f, ' % ra)  # NEW in WCS version
	full_data.write('%.4f, ' % dec)  # NEW in WCS version
	full_data.write('%.2f, ' % qso_x)
	full_data.write('%.2f' % qso_y)
	full_data.write('\n')

	if quiet:
		if foundz > 0:
			verified_sources.write(
				"%d \t%.1f \t%.1f \t%f  1 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\t%f\t%f\t%f\t%.1f\t%.1f\n" % (
					run_id, px, py, z, verr, used, gtemplate, found, ra, dec, verr, quality, z, foundz, foundm,
					get_impact(qso_x, qso_y, px, py, z), foundb))
			stat_good.write("%f %f %f\n" % (z, quality, peakratio))
		else:
			verified_sources.write("%d \t%.1f \t%.1f \t%f  1 \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\t-\t-\t-\t%.1f\n" % (
				run_id, px, py, z, verr, used, gtemplate, found, ra, dec, verr, quality,
				get_impact(qso_x, qso_y, px, py, z)))

		# verified_sources.write("%d \t%.1f \t%.1f \t%f %d \t%.1f \t%d \t%d \t%d\t%f\t%f\t%f\t%d\n" % (run_id,px,py,z,pcatalog[run_id],verr,used,gtemplate,found,ra,dec,verr,quality))
		verified_sources.flush()
		# if oiifound: ratios.write("%d %f\n" % (pcatalog[run_id],height2b/height2a))
		plt.close()
	del raw_flux
	del spec
	del data1
