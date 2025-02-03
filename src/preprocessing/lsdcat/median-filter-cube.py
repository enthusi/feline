#! /usr/bin/env python
#
# FILE:    median-filter-cube.py
# AUTHOR:  E. C. Herenz  / T. Urrutia
# DATE:    June 2012 / May 2020 / Feb 2022
# DESCR.:  Apply a median filter in spectral direction to the MUSE Datacube.
#          After that apply a small Gaussian filter in the spectral direction
#          to the MUSE Datacube and subtract the filtered data from the cube.
#          (do not alter the variances though, i.e. median and gauss filter
#          are only applied to the flux datacube)

import argparse
import multiprocessing
import os
# python standard library
import sys
import time
from datetime import datetime

# numpy / scipy / astropy
import numpy as np
from astropy.io import fits

# lsdcat/lib/line_em_funcs.py
import line_em_funcs as lef
import lsd_cat_lib
from median_filter_lib import medfilt_cube_para

#  macOS since Sierra uses "spawn"
try:
    multiprocessing.set_start_method('fork')
except RuntimeError:
    pass

__version__ = lsd_cat_lib.get_version()

# get string of the commandline that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command + ' ' + entity

starttime = time.time()
now = datetime.now()
get_timestring = lef.get_timestring

# parse command line
parser = argparse.ArgumentParser(
    description="""
    Subtract an in spectral direction
    median-filtered version of the datacube.
      Can be used to remove sources that have
    significant detectable continuum signal within the datacube.
    """)
# options visible to the user
parser.add_argument("fitscube",
                    type=str,
                    help="""
                    The FITS File that contains the flux- and
                     variances datacube
                    """)
parser.add_argument("-S", "--signalHDU", default='DATA', type=str,
                    help=""" HDU name or number (0-indexed) containing the
                     flux in fitscube
                    (default: 'DATA').  """)
parser.add_argument("-V", "--varHDU", default='STAT', type=str,
                    help=""" HDU name or number (0-indexed) containing
                     the variances in
                    fitscube (default: 'STAT').  Variances will just
                    be copied over to the output cube.  Use
                    --varHDU=-1 if this is not desired or if there is
                    no variance cube.  """)
parser.add_argument("-o", "--output", default=None,
                    help=""" Name of the output FITS file.  (default:
                    medfil_W+<width>+<fitscube>) """)
parser.add_argument("-W", "--width", default=151, type=int,
                    help=""" The median filter width in spectral pixels.
                      If an even number is
                    given, the actual filter value will be this number
                    +1. (Default: 151) """)
parser.add_argument("-G", "--gwidth", default=0, type=float,
                    help=""" The sigma of a Gaussian kernel for
                     in spectral pixels.  If not 0,
                    the median filtered spaxels will be smoothed with
                    a Gaussian kernel prior to subtraction.  This
                    prevents high-frequency leaks.  (Default: 0, i.e.
                    no Gaussian post filtering, recommended range: 5 -
                    10.) """)
parser.add_argument("--passes", default=1, type=int,
                    help="""Number of filtering iterations.
                      Subsequent iterations
                    (i.e. median filtering the median filtered data
                    again and again - 2 or 3 times) can sometimes lead
                    to better results (Default: 1).""")
parser.add_argument("--ao", action='store_true',
                    help=""" Fill notch-filter blocked region with
                     linear interpolated values.
                    This prevents an debiases the median filter
                    subtracted datacube close to the edges of the
                    notch filter.  See also the parameters
                    --notch_start, --notch_finish, and --notch_bound.
                    """)
parser.add_argument("--notch_start", type=float, default=5803.,
                    help=""" Start of the notch-filter blocked wavelength
                     range (in Angstrom,
                    default: 5803).  """)
parser.add_argument("--notch_end", type=float, default=5970.,
                    help=""" End of the notch-filter blocked wavelength
                     range (in Angstrom,
                    default: 5970).  """)
parser.add_argument("--notch_bound", type=int, default=30,
                    help=""" Size of window before- and after notch-filter
                     blocked region to
                    determine base points for the interpolation.  The
                    base point is calculated as the median over this
                    window.  """)
parser.add_argument("-t", "--num_cpu", default=None, type=int,
                    help=""" Number of CPU cores to be used
                     (default: all available).  """)
parser.add_argument("--savefiltered", action='store_true',
                    help=""" Also save the filtered data prior to subtraction.
                      Useful for
                    testing purposes.  """
                    )
parser.add_argument("--memmap", action='store_true',
                    help=""" Use memory mapping.  Reduces memory usage,
                     but might also reduce
                    speed of execution.  """)

args = parser.parse_args()

num_threads = args.num_cpu
width = args.width
gwidth = args.gwidth
fitscube = args.fitscube
output = args.output
memmap = True if args.memmap else False
num_pass = args.passes
signalHDU = lef.int_or_str(args.signalHDU)
# hdu indices are passed to the programm in 0-indexed
varHDU = lef.int_or_str(args.varHDU)  # notation (or HDU names)

# number of threads, either defined by user or use all CPUs
if num_threads is None:
    num_threads = multiprocessing.cpu_count()
# filter width must be uneven
if width % 2 == 0:
    width = width + 1
# set default output file name if not set
if output is None:
    output = 'medfil_W' + str(width) + '_' + os.path.basename(fitscube)

# read the input data / headers
print(fitscube + ': Reading input flux (HDU ' + str(signalHDU) + ') ... ' +
      get_timestring(starttime))
fitscube_hdu_list = fits.open(fitscube)
data, header = lef.read_hdu(fitscube, signalHDU, memmap=memmap)
primary_header = fitscube_hdu_list[0].header

# actual operations
if args.ao:
    xax = lef.wavel(header)
    notch_start_pix = np.argwhere(xax <= args.notch_start)[-1]
    notch_end_pix = np.argwhere(xax >= args.notch_end)[0]
    print(str(fitscube) + ': Linear interpolating over notch filter from ' +
          str(args.notch_start) + ' Angstrom (z_pix=' +
          str(notch_start_pix) + ') to ' +
          str(args.notch_end) + ' Angstrom (z_pix=' +
          str(notch_end_pix) + '.')
    print(str(fitscube) +
          ':  Width before & after filter to calculate median  ' +
          str(args.notch_bound) + ' (spectral bins).')
    pass_count = 0
    num_pass_total = num_pass
    data_in = data
    while num_pass > 0:
        pass_count += 1
        print(fitscube + ': Median filter - pass ' +
              str(pass_count) + '/' + str(num_pass_total))
        num_pass -= 1
        medfiltered_data = medfilt_cube_para(data_in, width, gwidth=gwidth,
                                             num_threads=num_threads,
                                             ao=True,
                                             notch_start_px=notch_start_pix,
                                             notch_end_px=notch_end_pix,
                                             notch_bound=args.notch_bound,
                                             verbose=True,
                                             fitscube=fitscube,
                                             starttime=starttime)
        data_in = medfiltered_data

    del data_in
else:
    pass_count = 0
    num_pass_total = num_pass
    data_in = data
    while num_pass > 0:
        pass_count += 1
        print(fitscube +
              ': Median filter - pass ' +
              str(pass_count) +
              '/' +
              str(num_pass_total))
        num_pass -= 1
        medfiltered_data = medfilt_cube_para(data_in, width, gwidth=gwidth,
                                             num_threads=num_threads,
                                             verbose=True, fitscube=fitscube,
                                             starttime=starttime)
        data_in = medfiltered_data

    del data_in

med_sub_data = data - medfiltered_data

# prepare for writing output to disk
primary_header = lef.mf_header(os.path.basename(sys.argv[0]),
                               __version__, width, gwidth, fitscube,
                               signalHDU, varHDU, primary_header)
primary_header['HISTORY'] = "Processed by LSDCat median-filter-cube.py " + \
                            " -- " + now.strftime("%m/%d/%Y, %H:%M:%S")
primary_header['HISTORY'] = '--- start of median-filter-cube.py command ---'
primary_header['HISTORY'] = command
primary_header['HISTORY'] = '--- end of median-filter-cube.py command ---'

if args.ao:
    primary_header['HIERARCH LSD MFRAO'] = (
        1, 'notch-gap interpolated')
    primary_header['HIERARCH LSD AO_S'] = (
        args.notch_start, '[Angstrom] notch gap start')
    primary_header['HIERARCH LSD AO_E'] = (
        args.notch_end, '[Angstrom] notch gap end')
    primary_header['HIERARCH LSD AO_W'] = (
        args.notch_bound, ' width of notch interp. (pixel)')

med_sub_header = header.copy()
med_sub_header['EXTNAME'] = \
    'MFGSUBDATA' if gwidth != 0 else 'MFSUBDATA'
med_header = header.copy()
med_header['EXTNAME'] = \
    'MFGDATA' if gwidth != 0 else 'MFDATA'
# write final output for object detection
primary_hdu = fits.PrimaryHDU(data=None, header=primary_header)
med_sub_data_hdu = fits.ImageHDU(
    data=med_sub_data.astype(np.float32), header=med_sub_header)

if varHDU != '-1':
    variance, varhead = lef.read_hdu(
        fitscube, varHDU, memmap=memmap)
    variance_hdu = fits.ImageHDU(
        data=variance.astype(np.float32), header=varhead)
    hdu_list = [primary_hdu, med_sub_data_hdu, variance_hdu]
else:
    hdu_list = [primary_hdu, med_sub_data_hdu]

if args.savefiltered:
    header['EXTNAME'] = 'MEDFIL'
    hdu_list.append(fits.ImageHDU(data=medfiltered_data.astype(np.float32),
                                  header=header))

# input cube contains exposure count HDU
try:
    ex_hdu = fitscube_hdu_list['EXP']
    hdu_list.append(ex_hdu)
except Exception as e:
    pass

out_hdu_list = fits.HDUList(hdus=hdu_list)
out_hdu_list.writeto(output, overwrite=True, output_verify='ignore')

print(str(fitscube) +
      ': All done! Filtered datacube'
      ' (with unchanged variances) stored in: ' +
      str(output) + ' (HDU0 - primary HDU header from input, ' +
      ' HDU 1: Original minus Medianf., (HDU2: Variances,'
      ' HDU3: Medianf. if desired) - ' +
      get_timestring(starttime))
