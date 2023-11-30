#! /usr/bin/env python
#
# FILE:   lsd_cc_spectral.py
# AUTHOR: Christian Herenz
# DESCR.: Second step of matched filtering: spectral cross-correlation.
# LICENSE: BSD 3-Clause License // https://opensource.org/licenses/BSD-3-Clause
#
# If you make use of this code in your research please cite:
# - Herenz, E. C., & Wisotzki, L. 2017,  A&A 602, A111.
#   https://doi.org/10.1051/0004-6361/201629507
# - Herenz, E. 2023, AN, e606
#   https://doi.org/10.1002/asna.20220091


import lib.lsd_cat_lib as lsd_cat_lib

__version__ = lsd_cat_lib.get_version()

import time
from datetime import datetime
import sys
import os
import numpy as np
from astropy.io import fits
import string
import random
import argparse
import multiprocessing
#  macOS since Sierra uses "spawn" 
multiprocessing.set_start_method('fork')

import gc

starttime = time.time()
now = datetime.now()

#########################################################################
# FUNCTIONS USED IN THIS PROGRAMM
    
from lib.wavelength_smooth_lib import *  # 

import lib.line_em_funcs as lef  # my own library with convenience functions
from lib.line_em_funcs import get_timestring
from lib.line_em_funcs import int_or_str

##########################################################################
# command line parsing

# get string of the commandline that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command+' '+entity

parser = argparse.ArgumentParser(description="""lsd_cc_spectral.py

Wavelength smoothing of all spaxels in the datacube with a gaussian
kernel. Operation is performed on signal and noise HDU.  """)
parser.add_argument("-i","--input",
                    type=str,
                    required=True,
                    help=""" Name of the input FITS file containing the spatially
                    cross-correlated datacube (output from lsd_cc_spatial.py).  """)
parser.add_argument("-v", "--varspec",
                    type=str,
                    help="""Name of the FITS file containing the variance spectrum.
                    (this option is ignored in --classic mode).
                    """ ,
                    default=None)
parser.add_argument("-e", "--expmap", type=str, default=None,
                    help=""" 
                    Name of the FITS file that contains the exposure map N_exp(x,y)
                    (default: none).  If specified, the matched-filter output
                    SN(x,y,z) will be weighted SN'(x,y,z) =
                    sqrt(N_exp(x,y)/max(N_exp)) · SN(x,y,z), where max(N_exp) will be
                    derived from the exposure map.  However, you can also fix this
                    value with the --nmax argument (usefull if it is not suitable to
                    extract the variance spectrum in the region of maximum depth).
                    Note that --exmpap is ignored if --classic mode is used.  """)
parser.add_argument("-F","--FWHM",
                    type=float,
                    default=300,
                    help="""Specify the FWHM of the Gaussian line template in km/s. [Default:
                    300 km/s] """)
parser.add_argument("-o","--output",
                    type=str,
                    default='',
                    help=""" Name of the output FITS file. The output FITS file will contain 2
                    HDUs: In HDU 0 the filtered signal is stored and
                    HDU 1 contains the propagated variances. [Default:
                    `cc_spec_+INPUT`,
                    i.e. `wavelength_smooth_` will be appended to the
                    input file name.  """)
parser.add_argument("-S",
                    "--SHDU",
                    type=int_or_str,
                    default='DATA_2DCC',
                    help=""" HDU name or number (0-indexed) in the input FITS file containing
                    the flux data. [Default: DATA_2DCC, i.e. the output HDU from lsd_cc_spatial.py] """)
parser.add_argument("-N",
                    "--NHDU",
                    type=int_or_str,
                    default=None,
                    help=""" In normal mode: HDU name or number in the variance spectrum FITS
                    file in which the 1D variances are stored; default
                    'VARSPEC'.  In classic mode: HDU name or number
                    (0-indexed) in the input FITS file containing the
                    variance data; default 'STAT_2DCC'.""")
parser.add_argument("-E", "--EHDU", type=int_or_str, default='EXP',
                    help="""HDU name or number (0-indexed) in the exposure map
                    FITS file where the exposure map is actually
                    stored.  (This parameter has no effect when
                    --classic mode is used).""")
parser.add_argument("-t",
                    "--threads",
                    type=int,
                    default=multiprocessing.cpu_count(),
                    help=""" Number of CPU cores used in parallel operation. [Default: all
                    available cpu cores] """)
parser.add_argument("--ignorenoise",
                    action="store_true",
                    help=""" Classic mode only: Switch to not propagate the variance.  If set
                    the output FITS file will contain only 1 HDU that stores the filtered
                    signal. (Note: In normal mode also one HDU is written to disk, but
                    this is correctly normalised to be interpreted as matched filter
                    output.) """)
parser.add_argument("--cunit3",
                    default='',
                    type=str,
                    help=""" Specify wavelength unit ('Angstrom' or 'nm' or 'um'). [Default:
                    Value from FITS Header.]  """)
parser.add_argument("--nanfile",
                    type=str,
                    default='none',
                    help=""" Name of an FITS file that contains a 2D image in` --nanmaskhdu`
                    (see below), that is of the same spatial dimensions as the input cube.
                    Spectra corresponding to NaNs in this image will be ignored in the
                    filtering. [Default: None] """)
parser.add_argument("--nanhdu",type=int_or_str,
                    default='4',
                    help=""" Number of HDU (0-indexed) or name of FITS file specified in
                    --namask, where the 2D image is stored. [Default: 4] """)
parser.add_argument("--classic", action='store_true',
                    help="""LSDCat 1.0 mode: Cross-correlate data with filter and propagate the
                    variance; cross-correlated data and propagated
                    varaince will be written to disk. This requires a variance cube. """)
parser.add_argument("--notemp", action='store_true',
                    help="""In --classic mode, do not write out temporary file to free memory
                    prior to propagating the variances.  This can speed up things if you have
                    a lot of memory.
""")
parser.add_argument("--nmax", type=float, default=0,
                    help=""" 
                    Weigh the resulting SN(x,y,z) cube via SN'(x,y,z) =
                    sqrt(N_exp(x,y)/nmax)) · SN(x,y,z), where N_exp(x,y) are the
                    values of the exposure map provided with the --expmap parameter.
                    If not set we use SN'(x,y,z) = sqrt(N_exp(x,y)/max(N_exp)) ·
                    SN(x,y,z) for the weighting (see --expmap).
""")
parser.add_argument("--overwrite", action='store_true',
                    help="Overwrite output files. Use with caution!")

# HIDDEN OPTION TO EXPERIMENT WITH MEMORY MAPPING
parser.add_argument("--memmap",action='store_true',help=argparse.SUPPRESS)


args = parser.parse_args()

inputfile = args.input
velocity = args.FWHM
num_threads = args.threads
classic = args.classic

if args.output == '':
    outfilename = 'cc_spec_' + inputfile
else:
    outfilename = args.output

data_hdu = args.SHDU     
if classic and args.NHDU == None:
    noise_hdu = 'STAT_2DCC'
elif not classic and args.NHDU == None:
    noise_hdu = 'VARSPEC'
else:
    noise_hdu = args.NHDU

# welcome message
program_name =  __file__.split('/')[-1]
print(program_name+' from LSDCat version '+str(__version__))
print(program_name+': Using ' +str(num_threads) + ' parallel threads')
if classic:
    print(program_name+': LSDCat classic mode!  Cross-correlatted data and propagated ' + \
          'variance will be written to disk.')

###########################################################################
# ACTUAL ROUTINE STARTS HERE

if args.memmap:
    print(inputfile+': Memory mapping in use (as requested with \
    --memmap)...')

print(inputfile+': Reading in datacube (HDU '+str(data_hdu)+')... '+\
      get_timestring(starttime))
# Nans in data are set to 0 (nans_to_value = True)
data, data_header = lef.read_hdu(inputfile,
                                 data_hdu,
                                 nans_to_value=True,
                                 memmap=args.memmap)

primary_header = fits.getheader(inputfile, 0)

if not classic and args.varspec != None:
    print(inputfile+': Reading input variance spectrum ' + args.varspec + \
          '(NHDU ' + str(noise_hdu) + ').')
    var_spec, var_header = fits.getdata(args.varspec, noise_hdu, header=True)
    if len(var_spec) != data.shape[0]:
        print('ERROR: Mismatch of spectral dimensions between variance spectrum (' + \
              str(len(var_spec)) + ') and input cube (' + str(data.shape[0]) + ').')
        sys.exit(2)
elif not classic and not args.ignorenoise and args.varspec == None:
    print('ERROR: No variance spectrum FITS file provided (-v), but --classic switch not set. ')
    sys.exit(2)
elif classic and args.varspec != None:
    print(inputfile+': --classic switch -- ignoring supplied 1D variance FITS file ' + \
          args.varspec + ' / 3D variance will be read from HDU '+str(noise_hdu)+' of ' + \
          inputfile + ' ...')
    
if args.cunit3 == '':
    # get values from header 
    cunit3 = data_header.get('CUNIT3')
else:
    cunit3 = args.cunit3
assert cunit3 == 'Angstrom' or cunit3 == 'nm' or cunit3 == 'um', \
    'unsupported wavelength unit!'

# 
if cunit3 == 'Angstrom':
    multiplier = 1
elif cunit3 == 'nm':
    multiplier = 10
elif cunit3 == 'um':
    multiplier = 1E4
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
startwav = (1 - crpix3) * cdelt3 +  crval3 
                  
# reshape the data for 1D-iteration 
print(inputfile+': Prepare the data for filtering... '+\
          get_timestring(starttime))
shape_0 = data.shape[0]
shape_1 = data.shape[1]
shape_2 = data.shape[2]
data = data.reshape(shape_0, shape_1 * shape_2)

# Spaxels to be ignored?
if args.nanfile != 'none':
    print(inputfile+': Using ' + str(args.nanfile) + ' (HDU:' + \
              str(args.nanhdu)+ ') to ignore NaN spaxels.') 
    print(get_timestring(starttime))
    nan_hdu = fits.open(args.nanfile)
    nans = nan_hdu[args.nanhdu].data
    assert nans.shape == (shape_1,shape_2)
    nans = nans.reshape(shape_1 * shape_2)
    nans_select = np.isnan(nans)
    num_nans = np.sum(nans_select)
    print(inputfile + ': Ignoring ' + str(num_nans) + \
          ' spaxels, because they are NaNs.')
else:
    print(inputfile+': No --nanfile set - all spaxels will be filtered... ' + \
              get_timestring(starttime))
    nans_select = None

print(inputfile+': Create the filter matrix... ' + \
          get_timestring(starttime))
# use the function create_filter_matrix from wavelength_smooth_lib.py
# to create the matrix that is used for cross correlation
# to search for emission line galaxies
if classic or args.ignorenoise:
    filter_matrix = create_filter_matrix_vel(velocity,
					     lambda_start=startwav,
					     cdelt=cdelt3,
					     lMax=shape_0)
    # square each entry of the matrix for error propgation
    filter_matrix_squared = np.square(filter_matrix)  
else:
    filter_matrix = create_filter_matrix_vel(velocity,
                                             lambda_start=startwav,
                                             cdelt=cdelt3,
                                             lMax=shape_0,
                                             var_spec=var_spec)

# filtering the flux
print(inputfile+': Threaded filtering starts... ' + \
          get_timestring(starttime))
filtered_data = filter_parallel(filter_matrix,
                                data,
                                num_threads,
                                filename=inputfile,
                                nans_select=nans_select)

# reshaping 2D list of spectra to 3D datacube again
filtered_data = filtered_data.reshape(filter_matrix.shape[0], shape_1, shape_2)
# truncation to actual spectral range:
start = int((filter_matrix.shape[0] - filter_matrix.shape[1])/2)
end = int(filtered_data.shape[0] - start)
filtered_data = filtered_data[start:end, :, :]

# primary header 
primary_header = lef.ccl_header(sys.argv[0], __version__, 
                                inputfile, data_hdu, noise_hdu, velocity,
                                primary_header, varspec=args.varspec)
primary_header['HISTORY'] = "Processed by LSDCat lsd_cc_spectral.py " + \
    " -- "+now.strftime("%m/%d/%Y, %H:%M:%S")
primary_header['HISTORY'] = '--- start of lsd_cc_spectral.py command ---'
primary_header['HISTORY'] = command
primary_header['HISTORY'] = '--- end of lsd_cc_spectral.py command ---'


# weight by exposure map, if desired
if args.expmap != None and not classic:
    print(inputfile + ': Weighing result with Nₑₓₚ(x,y) map from ' + \
          args.expmap + ' (HDU '+str(args.EHDU)+')')
    expmap, expmap_header = fits.getdata(args.expmap, args.EHDU, header=True)
    expmap_primary_header = fits.getheader(args.expmap, 0)
    assert (expmap.shape[0],
            expmap.shape[1]) == (filtered_data.shape[1],
                                 filtered_data.shape[2]), "expmap ↔ datacube shape missmatch"
    if args.nmax > 0:
        weights = np.sqrt(expmap / args.nmax)
    else:
        weights = np.sqrt(expmap / expmap.max())  
    filtered_data *= weights  # actual weighting

    # try to update header
    primary_header = lef.copy_hist_header(primary_header, expmap_primary_header)

if not classic and args.varspec != None:
    # update primary header with info from variance spectrum
    primary_header = lef.copy_hist_header(primary_header, var_header)

primary_hdu = fits.PrimaryHDU(data=None, header=primary_header)
data_header['EXTNAME'] = 'DATA_3DCC'
    
    
if not classic or args.ignorenoise:
    # done - write stuff to disk
    filtered_data_hdu = fits.ImageHDU(data=filtered_data, header=data_header)
    out_hdu_list = fits.HDUList(hdus=[primary_hdu, filtered_data_hdu])
    print(inputfile + ': Writing spectrally filtered data to disk: ' + \
          outfilename)

    try:
        out_hdu_list.writeto(outfilename, overwrite=args.overwrite, output_verify='ignore')
    except OSError:
        print('ERROR: Could not write output %s. Use --overwrite or different filename.' %outfilename)
        sys.exit(2)

    print(inputfile + \
          ': All done! Filtered datacube stored in : '+  str(outfilename) + ' ' + \
          get_timestring(starttime))

else:
    # Need to propagate the noise - same operation as above, but now
    # with squared filter matrix.  
    if not args.notemp:
        #  write out the filtered_data to a tempfile free some memory:
        tempfilename = ''.join(random.choice(string.ascii_uppercase + \
                                             string.digits) for x in range(6))
        tempfilename = tempfilename+'.fits'

        print(inputfile+': Writing temporary convolved flux datacube ('+\
              tempfilename+')... '+\
              get_timestring(starttime))
        lef.write_primary(filtered_data.astype(np.float32),
                          data_header,tempfilename)
    
        print(inputfile+': Freeing some memory... '+\
              get_timestring(starttime))
        del filtered_data
        del data
        gc.collect()
    
    print(inputfile+': Reading Noise (HDU '+str(noise_hdu)+')... '+\
              get_timestring(starttime))
    noise, noise_header = lef.read_hdu(inputfile,
                                       noise_hdu,
                                       nans_to_value=True,
                                       memmap=args.memmap)

    print(inputfile+': Prepare the noise for error propagation... '+\
              get_timestring(starttime))
    shape_1 = noise.shape[1]
    shape_2 = noise.shape[2]
    noise = noise.reshape(noise.shape[0], shape_1*shape_2)

    filtered_noise = filter_parallel(filter_matrix_squared,
                                     noise,
                                     num_threads,
                                     string="variance spectra",
                                     filename=inputfile,
                                     nans_select=nans_select)

    filtered_noise = filtered_noise.reshape(filter_matrix.shape[0],
                                            shape_1,shape_2)
    start = int((filter_matrix.shape[0] - filter_matrix.shape[1])/2)
    end = int(filtered_noise.shape[0] - start)
    filtered_noise = filtered_noise[start:end,:,:]

    noise_header['EXTNAME'] = 'STAT_3DCC'
   
    print(inputfile+': Preparing for writing out of the final datacube.... '+\
          get_timestring(starttime))

    if not args.notemp:
        convflux, convfluxhead = lef.read_hdu(tempfilename,0,
                                              memmap=args.memmap)
        os.remove(tempfilename)
    else:
        convflux, convfluxhead = (filtered_data.astype(np.float32),
                                  data_header)

    filtered_data_hdu = fits.ImageHDU(data=convflux.astype(np.float32),
                                      header=convfluxhead)
    filtered_noise_hdu = fits.ImageHDU(data=filtered_noise.astype(np.float32),
                                       header=noise_header)
    out_hdu_list = fits.HDUList(hdus=[primary_hdu,filtered_data_hdu,
                                      filtered_noise_hdu])


    print(inputfile+': Writing out the final datacube ... '+\
          get_timestring(starttime))
    out_hdu_list.writeto(outfilename, overwrite=args.overwrite, output_verify='ignore')

    print(inputfile + \
          ': All done! Wavelength smoothed cube & propagated noise stored in: '+\
          str(outfilename) + ' ' +  get_timestring(starttime))
