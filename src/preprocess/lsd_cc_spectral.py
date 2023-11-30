#! /usr/bin/env python
#
# FILE:   lsd_cc_spectral.py
# DESCR.: - convolve all spectra  with a gaussian filter with varying width          
#         - parallel (multiprocessing)
#         - with correct error propagation

__version__ = '1.0.2'

import time
import sys
import os
import math as m
import pylab as p
import numpy as np
from astropy.io import fits
import string
import random
import argparse
import multiprocessing
from scipy import signal

starttime = time.time()

#########################################################################
# FUNCTIONS USED IN THIS PROGRAMM

if sys.version_info[0] == 3:
    print('Incompatible to Python Versions > 3.0')
    sys.exit(2)


from lib.wavelength_smooth_lib import *  # my own libary

from lib.line_em_funcs import get_timestring 
import lib.line_em_funcs as lef # my own library with convenience functions

##########################################################################
# command line parsing

# get string of the commandline that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command+' '+entity

parser = argparse.ArgumentParser(description="""
lsd_cc_spectral.py

Wavelength smoothing of all spaxels in the datacube with a gaussian
kernel. Operation is performed on signal and noise HDU.
""")
parser.add_argument("-i","--input",type=str,required=True,help="""
Name of the input FITS file containing the flux (and variance) datacube.
""")
parser.add_argument("-F","--FWHM",type=float,default=300,help="""
Specify the FWHM of the Gaussian line template in km/s. [Default: 300 km/s]
""")
parser.add_argument("-o","--output",type=str,default='',help="""
Name of the output FITS file. The output FITS file will contain 2 HDUs: In HDU 0 the filtered signal is stored and HDU 1 contains the propagated variances. [Default: `wavelength_smooth_+INPUT`, i.e. `wavelength_smooth_` will be appended to the input file name.
""")
parser.add_argument("-S","--SHDU",type=int,default=0,help="""
HDU number (0-indexed) in the input FITS file containing the flux data. [Default: 0]
""")
parser.add_argument("-N","--NHDU",type=int,default=1,help="""
HDU number (0-indexed) in the input FITS file containing the variance data. [Default: 1]
""")
parser.add_argument("-t","--threads",type=int,
                    default=multiprocessing.cpu_count(),help="""
Number of CPU cores used in parallel operation. [Default: all available cpu cores]
""")
parser.add_argument("--ignorenoise",action="store_true",help="""
Switch to not propagate the variance.  If set the output FITS file will contain only 1 HDU that stores the filtered signal.
""")
parser.add_argument("--cunit3",default='',type=str,help="""
Specify wavelength unit ('Angstrom' or 'nm'). [Default: Value from FITS Header.]
""")
parser.add_argument("--memmap",action='store_true',help="""
Switch to use memory mapping.  Reduces memory usage, but also increases execution time.
""")
parser.add_argument("--nanfile",type=str,default='none',help="""
Name of an FITS file that contains a 2D image in` --nanmaskhdu` (see below), that is of the same spatial dimensions as the input cube.  Spectra corresponding to NaNs in this image will be ignored in the filtering. [Default: None]
""")
parser.add_argument("--nanhdu",type=int,default=4,help="""
Number of HDU (0-indexed) of FITS file specified in --namask, where the 2D image is stored. [Default: 4]
""")

args = parser.parse_args()

inputfile = args.input
velocity = args.FWHM
num_threads = args.threads

if args.output == '':
    outfilename = 'wavelength_smooth_'+inputfile
else:
    outfilename = args.output

memmap = True if args.memmap else False

data_hdu = args.SHDU 
noise_hdu = args.NHDU 

# TODO: welcome message (see lsd_cc_spatial.py)

###########################################################################
# ACTUAL ROUTINE STARTS HERE


if memmap:
    print(inputfile+': Memory mapping in use (as requested with \
    --memmap)...')

print(inputfile+': Reading Signal (HDU '+str(data_hdu)+')... '+\
      get_timestring(starttime))
# Nans in data are set to 0 (nans_to_value = True)
data,data_header = lef.read_hdu(inputfile,
                                data_hdu,
                                nans_to_value=True,
                                memmap=memmap)

# get values from header 
if args.cunit3 == '':
    cunit3 = data_header.get('CUNIT3')
else:
    cunit3 = args.cunit3
assert cunit3 == 'Angstrom' or cunit3 == 'nm'

# caclua
if cunit3 == 'Angstrom':
    multiplier = 1
elif cunit3 == 'nm':
    multiplier = 10
crval3 = data_header['CRVAL3']
try: 
    crpix3 = data_header['CRPIX3']
except KeyError:
    crpix3 = 1 # if no crpix3 value is given assume one
    print(inputfile+': No CRPIX3 value in header - assuming its 1!')
cdelt3 = data_header.get('CD3_3')
if not cdelt3:
    cdelt3 = data_header.get('CDELT3')
crval3 = crval3 * multiplier
cdelt3 = cdelt3 * multiplier
startwav = crval3 
                  
# reshape the data for 1D-iteration 
print(inputfile+': Prepare the data for filtering... '+\
          get_timestring(starttime))
shape_0 = data.shape[0]
shape_1 = data.shape[1]
shape_2 = data.shape[2]
data = data.reshape(data.shape[0],shape_1 * shape_2)

# Spaxels to be ignored?
if args.nanfile != 'none':
    print(inputfile+': Using '+args.nanfile+' (HDU:'+\
              str(args.nanhdu)+') to ignore NaN spaxels.')+\
              get_timestring(starttime)
    nan_hdu = fits.open(args.nanfile)
    nans = nan_hdu[args.nanhdu].data
    assert nans.shape == (shape_1,shape_2)
    nans = nans.reshape(shape_1 * shape_2)
    nans_select = p.isnan(nans)
    num_nans = p.sum(nans_select)
    print(inputfile+': Ignoring '+str(num_nans)+' spaxels, because they are NaNs.')
else:
    print(inputfile+': No --nanfile set - all spaxels will be filtered... '+\
              get_timestring(starttime))
    nans_select = None

print(inputfile+': Create the filter Matrix... '+\
          get_timestring(starttime))
# use the function create_filter_matrix from wavelength_smooth_lib.py
# to create the matrix that is used for cross correlation
# to search for Lyman Alpha emitters
filter_matrix = create_filter_matrix_vel(velocity,
					 lambda_start=startwav,
					 cdelt=cdelt3,
					 lMax=shape_0)
filter_matrix_squared = p.square(filter_matrix)  # <- for error propagation

# filtering the flux
print(inputfile+': Threaded filtering starts... '+\
          get_timestring(starttime))
filtered_data = filter_parallel(filter_matrix,
                                data,
                                num_threads,
                                filename=inputfile,
                                nans_select=nans_select)

# reshaping 2D list of spectra to 3D datacube again
filtered_data = filtered_data.reshape(filter_matrix.shape[0],shape_1,shape_2)
# truncation to actual spectral range:
start = (filter_matrix.shape[0]-filter_matrix.shape[1])/2
end = filtered_data.shape[0] - start
filtered_data = filtered_data[start:end,:,:]

# header
data_header['EXTNAME'] = 'FILTERED_DATA'
data_header['CCL'] = (sys.argv[0],
                      'spectral cross-correlation (CCL) routine')
data_header['CCLV'] = (__version__,
                       'CCL version')
data_header['CCLC'] = (command,
                       'CCL full command')
data_header['CCLIN'] = (inputfile,
                        'CCL input filename')
data_header['CCLINS'] = (data_hdu,
                         'CCL input data HDU - 0-indexed')
data_header['CCLINN'] = (noise_hdu,
                         'CCL input variance HDU - 0-indexed')
data_header['CCLVFWHM'] = (velocity,
                           'CCL filter FWHM [km/s]')


# now, if requested - same stuff as above is done for the noise
if not args.ignorenoise:
    # write out the filtered_data (temporary) to free some memory:
    tempfilename = ''.join(random.choice(string.ascii_uppercase +\
                                             string.digits) for x in range(6))
    tempfilename = tempfilename+'.fits'

    print(inputfile+': Writing temporary convolved flux datacube ('+\
          tempfilename+')... '+\
              get_timestring(starttime))
    lef.write_primary(filtered_data,data_header,tempfilename)
    print(inputfile+': Freeing some Memory... '+\
              get_timestring(starttime))
    del filtered_data
    print(inputfile+': Reading Noise (HDU '+str(noise_hdu)+')... '+\
              get_timestring(starttime))
    noise,noise_header = lef.read_hdu(inputfile,
                                      noise_hdu,
                                      nans_to_value=True,
                                      memmap=memmap)

    print(inputfile+': Prepare the noise for error propagation... '+\
              get_timestring(starttime))
    shape_1 = noise.shape[1]
    shape_2 = noise.shape[2]
    noise = noise.reshape(noise.shape[0],shape_1*shape_2)

    filtered_noise = filter_parallel(filter_matrix_squared,
                                     noise,
                                     num_threads,
                                     string="variance spectra",
                                     filename=inputfile,
                                     nans_select=nans_select)

    filtered_noise = filtered_noise.reshape(filter_matrix.shape[0],shape_1,shape_2)
    start = (filter_matrix.shape[0]-filter_matrix.shape[1])/2
    end = filtered_noise.shape[0] - start
    filtered_noise = filtered_noise[start:end,:,:]

    noise_header['EXTNAME'] = 'FILTERED_STAT'
    for key in data_header.iterkeys():
        if 'CCL' in key:
            noise_header[key] = data_header[key]
            
    convflux,convfluxhead = lef.read_hdu(tempfilename,0,
                                         memmap=memmap)
else:
    print(inputfile+': Ignoring noise on request - only filtered flux '+\
              'will be written to disk! '+\
              get_timestring(starttime))
    convflux,convfluxhead = filtered_data,data_header

print(inputfile+': Preparing for writing out of the final datacube.... '+\
          get_timestring(starttime))

filtered_data_hdu = fits.PrimaryHDU(data=convflux.astype(np.float32),
                                    header=convfluxhead)
if not args.ignorenoise:
    filtered_noise_hdu = fits.ImageHDU(data=filtered_noise.astype(np.float32),
                                       header=noise_header)
    out_hdu_list = fits.HDUList(hdus=[filtered_data_hdu,
                                      filtered_noise_hdu])
    os.remove(tempfilename)
else:
    out_hdu_list = fits.HDUList(hdus=[filtered_data_hdu])

out_hdu_list.writeto(str(outfilename),
                     clobber=True,output_verify='silentfix')

print(inputfile+\
          ': All done! Wavelength smoothed cube & propagated noise stored in: '+\
          str(outfilename)+' '+\
          get_timestring(starttime))
