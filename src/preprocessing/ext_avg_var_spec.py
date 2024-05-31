#! /usr/bin/env python
#
# FILE:   ext_avg_var_spec.py
# AUTHOR: Christian Herenz (2022)
# DESCR.: Extract average variance spectrum from variance HDU; needed as input
#         for lsd_cc_spectral.py

import argparse
import os
import sys
from datetime import datetime

import numpy as np
from astropy.io import fits

import line_em_funcs as lef  # my own library with convenience functions
import lsd_cat_lib
from line_em_funcs import int_or_str

__version__ = lsd_cat_lib.get_version()

# get string of the commandline that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command + ' ' + entity

now = datetime.now()

parser = argparse.ArgumentParser(description="""ext_avg_var_spec.py

Extract average variance spectrum from variance HDU; needed as input
for lsd_cc_spectral.py.

Make sure that the aperture covers a region where maximum number of
exposures, Nₘₐₓ, contributes to the datacube.  This is needed, as the
re-scaling with the 1D variance σ² to a variance map σ²(x,y) uses
σ²(x,y) = σ² · Nₘₐₓ) / Nₑₓₚ(x,y), where Nₑₓₚ(x,y) is the exposure
count map.  Obviously, this assumes that the integration time for each
exposure is the same.


""")
parser.add_argument("-i", "--input", type=str, required=True,
                    help="Name of input FITS file containin the variance HDU.")
parser.add_argument("-o", "--output", type=str,
                    help=""" Name of the output FITS file
                     containing the variance spectrum.
                    [Default: eavs_+INPUT, i.e. eavs_ will be
                    appended to the input file name]""")
parser.add_argument("-N", "--NHDU", type=int_or_str,
                    help="""HDU number (0-indexed) or name
                    in the input FITS file containing
                    the variance data. [Default: 'STAT']""", default='STAT')
parser.add_argument("-x", "--xcen", type=float,
                    help="""X cooordinate (1-indexed) of the central spaxel
                     of the extraction aperture.
                    [Default: Central spaxel along the x-axis of the
                    datacube]""", default=None)
parser.add_argument("-y", "--ycen", type=int,
                    help="""Y cooordinate (1-indexed) of the central spaxel
                     of the extraction aperture.
                    [Default: Central spaxel along the y-axis of the
                    datacube]""", default=None)
parser.add_argument("-r", "--radius", type=float,
                    help=""" Radius of the circular extraction aperture
                     (in spaxel).  [Default:
                    10] """, default=10.)
args = parser.parse_args()

if args.output is None:
    out_filename = 'eavs_' + args.input
else:
    out_filename = args.output

program_name = __file__.split('/')[-1]
print(program_name + ' from LSDCat version ' + str(__version__))
print(program_name + ' Extracting an average variance spectrum from ' +
      args.input +
      ' (HDU = ' + str(args.NHDU) + ')')

primary_header = fits.getheader(args.input, 0)
primary_header['HISTORY'] = "The average variance spectrum is from LSDCat2."
primary_header['HISTORY'] = ("ext_avg_var_spec.py -- " +
                             now.strftime("%m/%d/%Y, %H:%M:%S"))
primary_header['HISTORY'] = '--- start of ext_avg_var_spec.py command ---'
primary_header['HISTORY'] = command
primary_header['HISTORY'] = '--- end of ext_avg_var_spec.py command ---'

stat_cube, stat_header = lef.read_hdu(args.input, args.NHDU,
                                      nans_to_value=True, memmap=False)
# TODO: some spatial averaging here (maybe?!)

xcen = stat_cube.shape[2] / 2 if args.xcen is None else args.xcen - 1
ycen = stat_cube.shape[1] / 2 if args.ycen is None else args.ycen - 1

print(program_name + ': Extracting at (x_cen,y_cen) = (' +
      str(xcen) + ',' + str(ycen) + ')' +
      ' r = ' + str(args.radius) + ' px ...')
output_spec = lef.extract_spectra_circ(stat_cube, xcen, ycen, args.radius)
output_spec[output_spec == 0.] = np.nan

print(program_name + ' ... Done! Writing output to ' + out_filename)

# write out the avg variance spectrum
spec_hdu = fits.ImageHDU(data=output_spec)
spec_hdu.header['EXTNAME'] = 'VARSPEC'
spec_hdu.header['CRVAL1'] = stat_header['CRVAL3']
spec_hdu.header['CRPIX1'] = stat_header['CRPIX3']
spec_hdu.header['CDELT1'] = stat_header['CD3_3']
spec_hdu.header['CUNIT1'] = stat_header['CUNIT3']
spec_hdu.header['CTYPE1'] = stat_header['CTYPE3']
spec_hdu.header['BUNIT'] = stat_header['BUNIT']
spec_hdu.header['EAVSCMD'] = command

out_hdu_list = fits.HDUList([fits.PrimaryHDU(data=None, header=primary_header),
                             spec_hdu])
out_hdu_list.writeto(out_filename)
print(program_name + ' is done!')
