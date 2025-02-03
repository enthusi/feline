#! /usr/bin/env python

# FILE: s2n-cube.py
# AUTHOR: E. C. Herenz
# DESCRIPTION: Create a Signal-to-Noise Datacube from a datacube
#              with <Signal-HDU> and <Noise-HDU>.
#              SN = <Signal-HDU>/<Noise-HDU>
#              (if variances are stored in the the header, sqrt first!)
# (part of LSDCat)

import lsd_cat_lib

__version__ = lsd_cat_lib.get_version()

from line_em_funcs import int_or_str

import sys
import os
import argparse
from astropy.io import fits
import numpy as np
from datetime import datetime

import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)

now = datetime.now()

# store command that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command + ' ' + entity

parser = argparse.ArgumentParser(
    description="""Create a signal to noise datacube from a
     FITS file containing a signal and a noise HDU.""",
    epilog="(HDU numbering is 0-indexed)",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument("-i", "--input", required=True, help="""
Name of input FITS datacube file (mandatory argument).
""")
parser.add_argument("-n", "--nanmask", type=str, default=None, help="""
FITS file containing a 2D image of same spatial dimensions as cubes.
NaNs in this image are then a mask to exclude detections outside the
field of view. Best used with with whitelight image.
""")
parser.add_argument("--nanmaskhdu", type=int, default=4, help="""
HDU number (0-indexed) containing the NaN mask.
""")
parser.add_argument("-o",
                    "--output",
                    default=None,
                    help=""" Name of output S/N datacube
                     (default: s2n_+<input>). """)
parser.add_argument("-S", "--SHDU", type=int_or_str, default=0,
                    help=""" Name or number of HDU
                     (0-indexed) in input FITS datacube
                    containing the signal.""")
parser.add_argument("-N", "--NHDU", type=int_or_str, default=1,
                    help="""Name or number of HDU (0-indexed)
                     in input FITS datacube containing
                    the noise.""")
parser.add_argument("--sigma", action='store_true', help="""
Switch to interpret noise as sigma.
""")
parser.add_argument("--clobber", action='store_true', help="""
Overwrite already existing output (use at your own risk).
""")
parser.add_argument("--float64", action='store_true', help="""
Write out data in 64 bit.
""")

args = parser.parse_args()
inputfile = args.input
if args.output is None:
    outputfile = 's2n_' + inputfile
else:
    outputfile = args.output

signalHDU = args.SHDU
noiseHDU = args.NHDU
variance_switch = args.sigma
clobber = args.clobber
if args.nanmask is not None:
    nanmask_file = args.nanmask
    nanmask_hdu = args.nanmaskhdu

hdu = fits.open(inputfile)
primary_header = hdu[0].header

print(inputfile + ': Reading Signal (HDU:' +
      str(signalHDU) + ') ...')
signal2noise = hdu[signalHDU].data

print(inputfile + ': Reading Noise (HDU:' +
      str(noiseHDU) + ') ...')
noise = hdu[noiseHDU].data

if not variance_switch:
    print(inputfile + ': Calculating sigma from variance...')
    noise = np.sqrt(noise)

print(inputfile + ': Reading Header...')
header = hdu[signalHDU].header.copy()

print(inputfile + ': Dividing signal by noise...')
signal2noise /= noise

if args.nanmask is not None:
    nanmask = fits.getdata(nanmask_file, nanmask_hdu)
    assert nanmask.shape == signal2noise[0, :, :].shape
    print(inputfile + ': Applying Nan-Mask ' + nanmask_file + ' (HDU: ' +
          str(nanmask_hdu) + ')')
    nanmask[~np.isnan(nanmask)] = 1
    signal2noise *= nanmask

del noise

if not args.float64:
    print(inputfile + ': Converting result to 32-bit float...')
    signal2noise = signal2noise.astype(np.float32)

header['EXTNAME'] = 'SIGNALTONOISE'
header['BUNIT'] = '1'
print(inputfile + ': Writing ' + str(outputfile) + ' to disk...\n')

primary_header['HISTORY'] = "Processed by LSDCat s2n-cube.py " + \
                            " -- " + now.strftime("%m/%d/%Y, %H:%M:%S")
primary_header['HISTORY'] = '--- start of s2n-cube.py command ---'
primary_header['HISTORY'] = command
primary_header['HISTORY'] = '--- end of s2n-cube.py command ---'
hdu_list = fits.HDUList(hdus=[fits.PrimaryHDU(data=None,
                                              header=primary_header),
                              fits.ImageHDU(data=signal2noise,
                                            header=header)])
hdu_list.writeto(outputfile, overwrite=clobber)

del signal2noise
print(inputfile + ': Success!\n')
