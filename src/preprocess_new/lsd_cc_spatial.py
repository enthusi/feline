#! /usr/bin/env python
#
# FILE:   lsd_cc_spatial.py
# AUTHOR: Christian Herenz
# DESCR.: First step of matched filtering: spatial cross-correlation 
# LICENSE: BSD 3-Clause License // https://opensource.org/licenses/BSD-3-Clause
#
# If you make use of this code in your research please cite:
# - Herenz, E. C., & Wisotzki, L. 2017,  A&A 602, A111.
#   https://doi.org/10.1051/0004-6361/201629507
# - Herenz, E. 2023, AN, e606
#   https://doi.org/10.1002/asna.20220091


import lib.lsd_cat_lib as lsd_cat_lib
__version__ = lsd_cat_lib.get_version()

import sys,os,string,random,warnings,time
from datetime import datetime
import argparse
import multiprocessing
import math as m
import pylab as p
import numpy as np
from astropy.io import fits
from scipy import signal
import gc

#  macOS since Sierra uses "spawn" 
multiprocessing.set_start_method('fork')

starttime = time.time()
now = datetime.now()

######################################################################
# definitions of functions

# Memmap is now hard-coded False, as it leads to significant slow-down
# when being turned on
memmap = False

import lib.spatial_smooth_lib as spatial_smooth_lib  # spatial filtering functions
import lib.line_em_funcs as lef  # my own library with convenience functions
from lib.line_em_funcs import get_timestring
from lib.line_em_funcs import int_or_str

######################################################################
# command line parsing

# store command that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command+' '+entity

# workaround to allow for negative numbers not being treated as command-line arguments
# https://stackoverflow.com/a/21446783/2275260
for i, arg in enumerate(sys.argv):
  if (arg[0] == '-') and arg[1].isdigit(): sys.argv[i] = ' ' + arg

# selectors
outSel = False
selMask = False

parser = argparse.ArgumentParser(description="""
Spatial cross correlation of all the layers in the datacube with a
wavelength dependent FSF kernel (approximated by a Moffat or a Gaussian).""",
                                 epilog="""
Note:  
Polynomials for the PSF FWHM and the beta parameter are in the form 
p(lam) = sum_n p_n (lam - lam0)^n, where lambda is in Angstrom and lam0 
is specified via --lambda0.  For the PSF FWHM the p_n's are defined via 
-pc and are assumed to be in units of arcsec/(Angstrom)^n - where n is the
order of the coefficient.   For the dimensionless beta parameter the p_n's 
are defined alike via -bc in units of 1/(Angstrom)^n.
""")
parser.add_argument("-i","--input",
                    required=True,
                    type=str,
                    help="""
                    Name of the input FITS file containing the flux (and variance) datacube.
                    """)
parser.add_argument("-o","--output",
                    type=str,
                    default=None,
                    help="""
                    Name of the output FITS file. The output FITS file will contain 2
                    HDUs: In HDU 0 the filtered signal is stored and HDU 1 contains the
                    propagated variances. [Default: `cc_spat_+INPUT`,
                    i.e. `cc_spat_` will be appended to the input file name.]
                    """)
parser.add_argument("-S","--SHDU",
                    type=int_or_str,
                    default='1',
                    help="""
                    HDU name or number (0-indexed) in the input FITS file containing the
                    flux data. [Default: 1]
                    """)
parser.add_argument("-N","--NHDU",
                    type=int_or_str,
                    default='2',
                    help="""
                    HDU number (0-indexed) or name in the input FITS file containing the
                    variance data. [Default: 2]
                    """)
parser.add_argument("--std",
                    action='store_true',
                    help="""
                    Some noise cubes are std not variance (e.g. in KMOS).  If set the
                    input noise cube will be squared to make it variance.
                    """)
parser.add_argument("--ignorenoise",
                    action='store_true',
                    help="""
                    Switch to not propagate the variance.  If set the output FITS file
                    will contain only 1 HDU that stores the filtered signal.
                    """)
parser.add_argument("-m","--mask",
                    type=str,
                    default=None,
                    help="""
                    Name of a FITS file containing a mask. [Default: none]
                    """)
parser.add_argument("-M","--MHDU",
                    type=int_or_str,
                    default='0',
                    help="""
                    HDU number (0-indexed) or name of the mask within MASK-file.  The
                    mask is supposed to be a 2D array with the same spatial dimensions
                    containing only ones and zeros.  Spaxels corresponding to zero-valued
                    pixels in the mask will be set to zero prior and post the spatial
                    convolution operation. [Default: 1]
                    """)
parser.add_argument("-P","--pixscale",
                    type=float,
                    default=0.2,
                    help="""
                    Size of a spaxel in arcseconds. [Default: 0.2]
                    """)
parser.add_argument("-t","--threads",
                    type=int,
                    default=multiprocessing.cpu_count(),
                    help="""
                    Number of CPU cores used in parallel operation. 
                    [Default: all available cpu cores]
                    """)
parser.add_argument("-bc", default=['3.5'], nargs='*',
                    help="""
                    List of polynomial coefficients (whitespace seperated) for the beta parameter of
                    the Moffat profile, i.e. -bc b0 b1 b2 .. bn and for m>n bm == 0.
                    (default: 3.5).
                    """)
parser.add_argument('-pc', default=['0.8'], nargs='*',
                    help="""
                    List of polynomial coefficients (whitespace seperated) for the PSF FWHM,
                    i.e. -pc p0 p1 ... pn and for m>n pm == 0.
                    (default: 0.8)
                    """)
parser.add_argument("--lambda0",
                    type=float,
                    default=7050,
                    help="""
                    Zero-point of the polynomials, in Angstrom.  [Default: 7050 Angstrom]
                    """)
parser.add_argument("--gaussian",
                    action='store_true',
                    help="""
                    Switch to use a Gaussian profile instead of the default Moffat
                    profile as spatial filter profile. The -bc parameter(s) will be
                    ignored in that case.
                    """)
parser.add_argument("-T","--truncconstant",
                    default=8,
                    type=float,
                    help="""
                    Parameter controlling the truncation of the filter window: the filter
                    is truncated at T*WIDTH-PARAM - were WIDTH-PARAM = sigma for Gaussian
                    and FWHM for Moffat. [Default: 8]
                    """)
parser.add_argument("--classic", action='store_true',
                    help="""LSDCat 1.0 mode: Cross-correlate data with filter and propagate the
                    variance; cross-correlated data and propagated
                    varaince will be written to disk.  """)
parser.add_argument("--notemp", action='store_true',
                    help="""In --classic mode, do not write out temporary file to free memory
                    prior to propagating the variances.  This can speed up things if you have
                    a lot of memory.
""")
parser.add_argument("--overwrite", action='store_true',
                    help="Overwrite output files. Use with caution!")
args = parser.parse_args()

# fixed trunc_constant - defines where filter windows is truncated
# filter_window == 0 <=> r > trunc_constant*width_parameter + 1
# width_parameter = sigma for Gaussian - FWHM for Moffat
# - with this choice the Moffat filter is ~ twice as large as the Gaussian
# at same trunc-constant
trunc_constant = args.truncconstant
input_filename = args.input
data_hdu = args.SHDU
stat_hdu = args.NHDU
mask_hdu = args.MHDU
use_gaussian = args.gaussian
pix_scale = args.pixscale
num_threads = args.threads
std = args.std
ignorenoise = args.ignorenoise
classic = args.classic
if args.mask != None:
    selMask = True
    mask_filename = args.mask
else:
    mask_filename = 'no_mask'
if args.output == None:
    out_filename = 'cc_spat_'+input_filename
else:
    out_filename = args.output


# polynomial coefficients for FWHM and beta
# (reversed, for np.polyval)
pc = [float(pcn) for pcn in args.pc][::-1]   
bc = [float(bcn) for bcn in args.bc][::-1]  
lambda_0 = args.lambda0

######################################################################
# program starts here

# welcome message:
program_name =  __file__.split('/')[-1]
filter_name = 'Gaussian' if use_gaussian else 'Moffat'
print(program_name+' from LSDCat version '+str(__version__))
print(program_name+' run on datacubes in inputfile: '+input_filename+\
          ' (flux in HDU: '+str(data_hdu)+', variance in HDU: '+\
          str(stat_hdu)+')')
if selMask == False:
    print(program_name+': Using no mask cube!')
else:
    print(program_name+': Using mask cube in '+mask_filename+' (HDU: '+\
              str(mask_hdu)+')')
print(program_name+': PSF shape model = '+filter_name)
print(program_name+': PSF FHWM(lambda) dependence via polynomial with coefficients ' + \
      ', '.join(['p'+str(i[0])+'='+str(i[1]) for i in enumerate(pc[::-1])]))
if filter_name == 'Moffat':
    print(program_name+': Moffat beta(lambda) dependence via polynomial with coefficients ' + \
      ', '.join(['b'+str(i[0])+'='+str(i[1]) for i in enumerate(bc[::-1])]))
print(program_name+': Reference wavelength lambda0 in PSF polynomial(s): ' + str(lambda_0) + \
      ' Angstrom')
print(program_name+': plate scale = '+str(pix_scale)+' arcsec')
print(program_name+': Using ' +str(num_threads) + ' parallel threads')
print(program_name+': Spatial covolution using Fast Fourier transform!')
if classic:
    print(program_name+': LSDCat classic mode!  Cross-correlatted data and propagated ' + \
          'variance will be written to disk.')

# reading in the data & headers
primary_header = fits.getheader(input_filename, 0)

print(input_filename+': Reading in the Data Cube... (HDU '+str(data_hdu)+') '+\
          get_timestring(starttime))
cube_data, cube_header = lef.read_hdu(input_filename,
                                      data_hdu,
                                      nans_to_value=False,
                                      memmap=memmap)  #first  we keep the nans, but
# we set all NaNs to zero here, otherwise the FFT filtering produces rubbish
nancube = np.isnan(cube_data)   # this nan-selector cube is also used
                                # in the end to set all original nans
                                # to zero (so FFT artifacts in the border regions are ignored)
cube_data[nancube] = 0  # note: with newer astropy > 3, convolve(...,nan_treatment='fill')
                        # might do the same (TODO: check)

crval3 = cube_header['CRVAL3']
try:
    cdelt3 = cube_header['CD3_3']
except KeyError:
    cdelt3 = cube_header['CDELT3']
crpix3 = cube_header['CRPIX3']
try:
    cunit3 = cube_header['CUNIT3']
except KeyError:
    cunit3 = 'Angstrom'
    print('WARNING: No CUNIT3 Keyword specifying the wavelength unit found. Assuming Angstrom.')

if selMask == True:
    print(str(input_filename)+': Applying mask '+str(mask_filename)+\
              ' (HDU'+str(mask_hdu)+') to Data (File '+str(input_filename)+\
              ', HDU '+str(data_hdu)+')... '+get_timestring(starttime))
    mask,mask_head = lef.read_hdu(mask_filename,
                                  mask_hdu,
                                  memmap=memmap)
    cube_data *= mask
else:
    mask = [] # just a dummy
    print('No mask is given... Performing all operations on umasked cube.\
 (File '+input_filename+', HDU '+str(data_hdu)+')... '+\
              get_timestring(starttime))

print(input_filename+': Threaded Filtering starts...'+\
          get_timestring(starttime))
print(input_filename+' ... Filter window PSFs are '+\
          filter_name+'s ...') # hope this is always the right plural ;-)

##################### ACTUAL COMPUTATION STARTS HERE ####################

# convolving all sclices (split over the processors) of data 
# with wavelength dependent gaussian kernel:

# 1. MAKE WAVELENGTH DEPENDENT FILTER PROFILES
# --------------------------------------------
length = cube_data.shape[0]
print(input_filename+\
          ': Creating the wavelength dependent PSF filter for '+\
          str(length)+' datacube layers. '+\
          get_timestring(starttime))

# xax is the array containing the wavelength values in Angstrom! the
# unit is important here - since a polynomial can be used to describe
# the wavelength dependence of the PSFs FWHM - in this case the
# coefficents describe a polynomial FWHM(lambda[Angstrom])
xax = crval3 + cdelt3 * (np.arange(length) + 1 - crpix3)  # Angstrom - default in MUSE
if cunit3 == 'nm':
    xax *= 10  # e.g. early MUSE cubes
elif cunit3 == 'um':
    xax *= 1E4  #  KMOS pipeline default
elif cunit3 != 'Angstrom':
    print('ERROR: Unknown wavelength unit: '+cunit3)
    sys.exit(2)

if use_gaussian:
    filter_windows = \
        spatial_smooth_lib.makeGaussians_poly(
            xax=(xax - lambda_0), p=pc, pix_scale=pix_scale,
            trunc_constant=trunc_constant, classic=classic)
else:
    filter_windows = \
        spatial_smooth_lib.makeMoffats_poly(
            xax=(xax - lambda_0), p=pc, b=bc, pix_scale=pix_scale,
            trunc_constant=trunc_constant, classic=classic)

print(input_filename + \
          ": Average size of the filter windows  " + \
          str(filter_windows[int(len(filter_windows)/2)].shape[0]) + \
          "^2 px. " + \
          get_timestring(starttime))

# create also the filters for the error propagation (below)
filter_windows_squared = [np.square(filter_window)
                          for filter_window in filter_windows]

# 2. ITERATATE OVER THE SLICES AND FILTER THEM 
# --------------------------------------------
# (split up over several processors)
filtered = spatial_smooth_lib.filter_parallel(filter_windows,
                                              cube_data,
                                              num_threads,
                                              selMask,
                                              mask,
                                              filename=input_filename,
                                              method='fft')
filtered[nancube] = 0  # setting all nans in the original cube to zero
                       # in the output datacube

mask = [] # delete old mask, and dummy it

# add information about processing with this script to the primary header
primary_header = lef.ccs_header(os.path.basename(sys.argv[0]),
                                __version__, input_filename,
                                data_hdu, stat_hdu, pix_scale,
                                mask_filename, mask_hdu, filter_name,
                                trunc_constant, pc, lambda_0, bc,
                                primary_header)
primary_header['HISTORY'] = "Processed by LSDCat lsd_cc_spatial.py " + \
    " -- "+now.strftime("%m/%d/%Y, %H:%M:%S")
primary_header['HISTORY'] = '--- start of lsd_cc_spatial.py command ---'
primary_header['HISTORY'] = command
primary_header['HISTORY'] = '--- end of lsd_cc_spatial.py command ---'

cube_header['EXTNAME'] = 'DATA_2DCC'

if not classic or ignorenoise:
    # we're done - write stuff to disk
    print(input_filename+': Writing out spatially cross-correlated flux ... '+\
          get_timestring(starttime))
    lef.write_fitscube(filtered, cube_header, primary_header, out_filename,
                       overwrite=args.overwrite)
    print(input_filename + ': All Done! Written spatial cross-corelated flux to '+\
          out_filename+' . '+ get_timestring(starttime))

else:
    # we need to propagte the variances .... (classic LSDCat)
    if not args.notemp:
        tempfilename = ''.join(random.choice(string.ascii_uppercase + string.digits) 
                               for x in range(6)) 
        tempfilename = tempfilename + '.fits'
        print(input_filename+': Writing out temporary filtered cube ' + \
              tempfilename+' ... ' + get_timestring(starttime))
        # writing out temporary FITS file, since we do not want to store 2
        # cubes (data and variance) in memory
        lef.write_primary(filtered, cube_header, tempfilename)
        del filtered
        del cube_data
        gc.collect()

    print(input_filename+': Reading in the Noise Cube... (HDU '+str(stat_hdu)+') '+\
              get_timestring(starttime))
    stat_data,stat_head = lef.read_hdu(input_filename,
                                       stat_hdu,
                                       nans_to_value=True,
                                       memmap=memmap)

    if std:
        # If the input noise cube is std (e.g., KMOS default), we need to square it
        stat_data = stat_data ** 2.

    stat_head['EXTNAME'] = 'STAT_2DCC'
    print(input_filename+': Error propagation in the noise cube... '+\
              get_timestring(starttime))

    # error propagation on all sclices of stat for the filtering
    filtered_stat = spatial_smooth_lib.filter_parallel(filter_windows_squared,
                                                       stat_data,
                                                       num_threads,
                                                       selMask,
                                                       mask,
                                                       filename=input_filename,
                                                       method='fft')

    print(input_filename + ': Writing out final data... ' + \
          get_timestring(starttime))
    if not args.notemp:
        filtered, filtered_header = lef.read_hdu(tempfilename,0,memmap=memmap)
        os.remove(tempfilename)
    else:
        filtered, filtered_header = filtered, cube_header
    filtered_data_hdu = fits.ImageHDU(data=filtered, header=filtered_header)

    del filtered
    del filtered_header
    filtered_noise_hdu = fits.ImageHDU(data=filtered_stat, header=stat_head)
    out_hdu_list = fits.HDUList(hdus=[fits.PrimaryHDU(data=None,
                                                      header=primary_header),
                                      filtered_data_hdu, filtered_noise_hdu])
    out_hdu_list.writeto(out_filename,
                         overwrite=args.overwrite,
                         output_verify='ignore')

print(input_filename + ': All Done: lsd_cc_spatial.py output written to ' + \
      out_filename+' . ' + get_timestring(starttime))
