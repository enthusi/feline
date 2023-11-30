#! /usr/bin/env python
#
# FILE:   lsd_cc_spatial.py
# AUTHOR: Christian Herenz (2011,2012,2013,2014,2015)
# DESCR.: Convolve all spectral slices with a gaussian filter, 
#         with FWHM given in arcseconds
#         with correct error propagation
#         uses pythons multiprocessing to split the work over 
#         several cores in a multicore machine
# CHANGELOG:
# 1.0alpha
# - intial version
# 1.0.1
# - PSF lambda dependence: coeffients in [''/Angstrom^n]
#   and p(lambda[A]) = p_0 + SUM_n p_n (lambda - lambda_0)^n

__version__ = '1.0.2'

import sys,os,string,random,warnings,time
import argparse
import multiprocessing
import math as m
import pylab as p
import numpy as np
from astropy.io import fits
from scipy import signal
from scipy import version

scipy_version = version.version
starttime = time.time()
 
if sys.version_info[0] == 3:
    print('Incompatible to Python Versions > 3.0')
    sys.exit(2)

######################################################################
# definitions of functions

assert scipy_version[1] >= 9 or scipy_version[0] >= 1

# Python versions < 2.7 memmap conflicts in multiprocessing fix:
memmap = False if sys.version_info[1] < 7 else True

import lib.line_em_funcs as lef # my own library with convenisky_blocks=\nce functions
import lib.spatial_smooth_lib as spatial_smooth_lib # <- library for
                                                    # making the
                                                    # spatial
                                                    # smoothing
                                                    # filtera with
                                                    # function
                                                    # makeGaussians
from lib.line_em_funcs import get_timestring 

######################################################################
# command line parsing

# get string of the commandline that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command+' '+entity

# convolve2d returns warning on complex recasting, 
# working only with real values - irritating thus disabled:
# warnings.simplefilter("ignore", np.ComplexWarning) 

# selectors
outSel = False
selMask = False

parser = argparse.ArgumentParser(description="""
Spatial cross correlation of all the layers in the datacube with a
wavelength dependent FSF kernel (approximated by a Moffat - can be
changed to Gaussian if requested) kernel. Operation is performed on 
signal and error is propagated accordingly.""",
                                 epilog="""
Note: The PSF widths refer always to the FWHM of the profile! 
The polynomial coefficients are assumed to be in units of 
FWHM/(Angstrom)^n - where n is the order of the coefficient.
""")
parser.add_argument("-i","--input",required=True,type=str,
                    help="Name of the input FITS file containing the flux (and variance) datacube.")
parser.add_argument("-o","--output",type=str,default=None,
                    help="Name of the output FITS file. The output FITS file will contain 2 HDUs: In HDU 0 the filtered signal is stored and HDU 1 contains the propagated variances. [Default: `spatial_smoothed_+INPUT`, i.e. `spatial_smoothed_` will be appended to the input file name.]")
parser.add_argument("-S","--SHDU",type=int,default=0,
                    help="HDU number (0-indexed) in the input FITS file containing the flux data. [Default: 0]")
parser.add_argument("-N","--NHDU",type=int,default=1,
                    help="HDU number (0-indexed) in the input FITS file containing the variance data. [Default: 1]")
parser.add_argument("--ignorenoise",action='store_true',
                    help="Switch to not propagate the variance.  If set the output FITS file will contain only 1 HDU that stores the filtered signal.")
parser.add_argument("-m","--mask",type=str,default=None,
                    help="Name of a FITS file containing a mask. [Default: none]")
parser.add_argument("-M","--MHDU",type=int,default=0,
                    help="HDU number (0-indexed) of the mask within MASK-file.  The mask is supposed to be a 2D array with the same spatial dimensions containing only ones and zeros.  Spaxels corresponding to zero-valued pixels in the mask will be set to zero prior and post the spatial convolution operation. [Default: 1]")
parser.add_argument("-P","--pixscale",type=float,default=0.2,
                    help="Size of a spaxel in arcseconds. [Default: 0.2]")
parser.add_argument("-t","--threads",type=int,default=multiprocessing.cpu_count(),
                    help="Number of CPU cores used in parallel operation. [Default: all available cpu cores]")
# parser.add_argument("--nofft",action='store_true',help="Don't use FFT convolution. Much slower!")
parser.add_argument("-b","--beta",type=float,default=3.5,
                    help="BETA parameter of the Moffat profile.  [Default: 3.5].")
parser.add_argument("--gaussian",action='store_true',
                    help="Switch to use a Gaussian profile instead of the default Moffat profile as spatial filter profile. The &beta; parameter will be ignored in that case.")
parser.add_argument("-p0",type=float,default=0.8,
                    help="0th order coefficient (in arcseconds) for polynomial approximation for PSF FWHM-lambda dependency. [Default: 0.8]")
parser.add_argument("-p1",type=float,default=0.0,
                    help="1st order polynomial coefficient (in arcseconds/Angstrom). [Default: 0.8]")
parser.add_argument("-p2",type=float,default=0.0,
                    help="2nd order polynomial coefficient (in arcsec/Angstrom**2). [Default: 0]")
parser.add_argument("--lambda0",type=float,default=7050,
                    help="Zero-point of the polynomial, in Angstrom.  [Default: 7050 Angstrom]")
parser.add_argument("-T","--truncconstant",default=8,type=float,
                    help="Parameter controlling the truncation of the filter window: the filter is truncated at T*WIDTH-PARAM - were WIDTH-PARAM = sigma for Gaussian and FWHM for Moffat. [Default: 8]")

args = parser.parse_args()

# if not args.nofft:
#     method = 'fft'
# else:
#     method = 'normal'

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

if args.mask != None:
    selMask = True
    mask_filename = args.mask
else:
    mask_filename = 'no_mask'
if args.output == None:
    out_filename = 'spatial_smoothed_'+input_filename
else:
    out_filename = args.output

#fwhm = args.FWHM 
beta = args.beta # beta is fixed
use_gaussian = args.gaussian
#air_mass = args.airmass
#scale_length = args.scalelength

pix_scale = args.pixscale
num_threads = args.threads

ignorenoise = args.ignorenoise

pc = [args.p0,args.p1,args.p2][::-1] # <- i.e. reversed, for np.polyval
lambda_0 = args.lambda0

######################################################################
# program starts here

# welcome message:
program_name =  __file__.split('/')[-1]
filter_name = 'Gaussian' if use_gaussian else 'Moffat'
print(program_name+' version '+__version__)
print(program_name+' run on datacubes in inputfile: '+input_filename+\
          ' (flux in HDU: '+str(data_hdu)+', variance in HDU: '+\
          str(stat_hdu)+')')
if selMask == False:
    print(program_name+': Using no mask cube!')
else:
    print(program_name+': Using mask cube in '+mask_filename+' (HDU: '+\
              str(mask_hdu)+')')
print(program_name+': PSF shape model = '+filter_name)
if filter_name == 'Moffat':
    print(program_name+': Moffat beta = '+str(beta))
psf_lambda_string = \
    'polynomial approximation using the coefficients [p_0,p_1,p_2] '+\
    str(pc[::-1])

print(program_name+': PSF lambda dependence via '+psf_lambda_string)
print(program_name+': Zero-Wavelength in Polynomial: '+str(lambda_0)+' Angstrom')
print(program_name+': Using '+str(num_threads)+' parallel threads')
print(program_name+': Spatial covolution using Fast Fourier transform!')

# reading in the flux data & header
# sets NaN values to 0, if nans_to_value = True
print(input_filename+': Reading in the Data Cube... (HDU'+str(data_hdu)+') '+\
          get_timestring(starttime))
cube_data,cube_header = lef.read_hdu(input_filename,
                                     data_hdu,
                                     nans_to_value=False,
                                     memmap=memmap)

# we set all NaNs to zero here, otherwise the FFT filtering produces rubbish
nancube = np.isnan(cube_data)   # this nan-selector cube is also used
                                # in the end to set all original nans
                                # to zero (so FFT artifacts in the border regions are ignored)
cube_data[nancube] = 0


print cube_data.shape
"""
print "* Martin Wendt: blanking sky planes now!... *"
sky_blocks=[\
  5577.3500,\
6300.3100,\
8827.2504,\
8827.1937,\
8919.5433,\
8919.7480,\
8885.8389,\
8885.9440,\
8344.6849,\
8344.6287,\
8957.9129,\
8958.2390,\
8430.0165,\
8430.2110,\
6363.7800,\
8399.0908,\
8399.1920,\
8465.1227,\
8465.4287,\
8778.5398,\
8778.4913,\
8836.6420,\
8836.4323,\
9001.0963,\
7913.6592,\
7913.6037,\
9001.5601,\
8768.1057,\
8768.0408,\
5889.9600,\
7993.0458,\
7993.2313,\
7750.6498,\
7750.6131,\
8903.3333,\
8903.3925,\
7821.3361,\
7821.4923,\
7964.4405,\
7964.5381,\
8943.6836,\
8943.6727,\
8299.0287,\
8298.9789,\
7276.4434,\
7276.4069,\
8352.9953,\
8352.7886,\
7853.1305,\
7853.3880,\
7794.0066,\
7794.0836,\
7340.7473,\
7340.8934,\
8504.5269,\
8504.9575,\
8867.7319,\
8867.8325,\
8025.4488,\
8025.7375,\
8288.7185,\
8288.6497]

lmin=4749.794921875
dpix=4
step=1.25

for pos in sky_blocks:
  
  pix=int((pos-lmin)/step+0.5)
  start=pix-dpix
  end=pix+dpix
  print pos,pix
  #print start,end,end-start
  #set sky-area to significance in S/N cube!
  #this way a line in this area will neither be good nor bad
  cube_data[start:end,:,:]=0.0
print "* Martin Wendt: blanking sky planes done... *"
"""
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

# 1. MAKE WAVELENGTH DEPENDENT FILTER SHAPES
# ------------------------------------------
length = cube_data.shape[0]
print(input_filename+\
          ': Creating the wavelength dependent PSF filter for '+\
          str(length)+' datacube layers. '+\
          get_timestring(starttime))

# xax is the array containing the wavelength values in Angstrom! the
# unit is important here - since a polynomial can be used to describe
# the wavelength dependence of the PSFs FWHM - in this case the
# coefficents describe a polynomial FWHM(lambda[Angstrom])
if cunit3 == 'nm':
    xax = 10*(crval3+cdelt3*(p.arange(length)-crpix3))
elif cunit3 == 'Angstrom':
    xax = crval3+cdelt3*(p.arange(length)-crpix3)
else:
    print('ERROR: Unknown wavelength unit: '+cunit3)
    sys.exit(2)

print(input_filename+\
          ": Using polynomial fit to approximate wavelength"+\
          " dependence of PSF FWHM ...")
print(input_filename+\
          ": Coefficients of the polynomial [p_0,p_1,p_2]="+str(pc[::-1])+\
          " (plate scale: "+str(pix_scale)+" arcsec^2 per spaxel)")
if use_gaussian:
    filter_windows_poly = spatial_smooth_lib.makeGaussians_poly
    filter_windows = filter_windows_poly(xax - lambda_0,
                                         pc,
                                         pix_scale=pix_scale,
                                         trunc_constant=trunc_constant)
else:
    filter_windows_poly = spatial_smooth_lib.makeMoffats_poly
    filter_windows = filter_windows_poly(xax - lambda_0,
                                         pc,beta,
                                         pix_scale=pix_scale,
                                         trunc_constant=trunc_constant)

print(input_filename+\
          ": Average size of the filter windows  "+\
          str(filter_windows[len(filter_windows)/2].shape[0])+\
          "^2 px. "+\
          get_timestring(starttime))

# create also the filters for the error propagation (below)
filter_windows_squared = []
for filter_window in filter_windows:
    filter_windows_squared.append(np.square(filter_window))

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
filtered[nancube] = 0  # setting all nans in the original cube to zero in the output datacube 

mask = [] # delete old mask, and dummy it

# header
cube_header['EXTNAME'] = 'DATA_2DSMOOTH'
cube_header['CCS'] = (os.path.basename(sys.argv[0]),
                      'spatial cross-correlation (CCS) routine')
cube_header['CCSV'] = (__version__,
                       'CCS version')
cube_header['CCSC'] = (command,
                       'CCS full command')
# cube_header['CCSM'] = (method,
#                        'CCS method - normal or fft')
cube_header['CCSIN'] = (input_filename,
                        'CCS input filename')
cube_header['CCSINS'] = (data_hdu,
                         'CCS input data HDU - 0-indexed')
cube_header['CCSINN'] = (stat_hdu,
                         'CCS input variance HDU - 0-indexed')
cube_header['CCSPXSC'] = (pix_scale,
                          'assumed pix scale arcsec/pix in datacube')
cube_header['CCSMSK'] = (mask_filename,
                         'SSC mask filename')
cube_header['CCSMSKH'] = (mask_hdu,
                          'CCS mask filename HDU')
cube_header['CCSFILT'] = (filter_name,
                          'CCS filter funtion - i.e. Moffat or Gaussian')
cube_header['CCSFILTT'] = (trunc_constant,
                           'CCS filter truncated at CCSFILTT x FWHM')
if filter_name == 'Moffat':
    cube_header['CCSMFB'] = (beta,
                             'CCS Moffat beta')
cube_header['CCSPLY'] = (True,
                         'poly p0+p1*(l-l0)+p2*(l-l0)**2')
cube_header['CCSPLYP0'] = (args.p0,
                           'p0 [arcec]')
cube_header['CCSPLYP1'] = (args.p1,
                           'p1 [arcsec/AA]')
cube_header['CCSPLYP2'] = (args.p2,
                           'p2 [arcsec/AA**2]')
cube_header['CCSPLYL0'] = (lambda_0,
                           'l0')

if not ignorenoise:
    # writing out temporary FITS file, since we do not want to store 2
    # cubes (data and variance) in memory
    # genrating random 6 digit or letter input_filename for temporary
    # fits file (e.g.  KNX9RC):
    tempfilename = ''.join(random.choice(string.ascii_uppercase + string.digits) 
                           for x in range(6))
    tempfilename = tempfilename+'.fits'
    print(input_filename+': Writing out temporary filtered cube '+tempfilename+' ... '+\
              get_timestring(starttime))
    lef.write_primary(filtered,cube_header,tempfilename)
    del filtered

    # NOW DOING THE SAME FOR THE NOISE, BUT WITH SQUARED MATRIX & WITHOUT MASK
    # (assuming noise is stored as variance)
    print(input_filename+': Reading in the Noise Cube... (HDU'+str(stat_hdu)+') '+\
              get_timestring(starttime))
    stat_data,stat_head = lef.read_hdu(input_filename,
                                       stat_hdu,
                                       nans_to_value=True,
                                       memmap=memmap)

    stat_head['EXTNAME'] = 'STAT_2DSMOOTH'
    for key in cube_header.iterkeys():
        if 'CCS' in key:
            stat_head[key] = cube_header[key]
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

    print(input_filename+': Writing out final data & deletion of temporary files... '+\
              get_timestring(starttime))
    filtered,filtered_header = lef.read_hdu(tempfilename,0,memmap=memmap)
    filtered_data_hdu = fits.PrimaryHDU(data=filtered,header=filtered_header)
    os.remove(tempfilename)
    del filtered ; del filtered_header
    filtered_noise_hdu = fits.ImageHDU(data=filtered_stat,header=stat_head)
    out_hdu_list = fits.HDUList(hdus=[filtered_data_hdu,filtered_noise_hdu])
    out_hdu_list.writeto(out_filename,clobber=True,output_verify='silentfix')
    print(input_filename+\
          ': All Done! Written spatial smoothed data & propagated error to '+\
          out_filename+' . '+\
          get_timestring(starttime))
else:
    print(input_filename+': Writing out cross-correlated flux... '+\
              get_timestring(starttime))
    cube_header['EXTNAME'] = 'DATA_2DSMOOTH'
    lef.write_primary(filtered,cube_header,out_filename)
    print(input_filename+\
          ': All Done! Written spatial smoothed data '+\
              out_filename+' . '+\
              get_timestring(starttime))
    
