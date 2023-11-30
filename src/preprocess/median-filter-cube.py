#! /usr/bin/env python
#
# FILE:    median-filter-cube.py
# AUTHOR:  C. Herenz
# DATE:    June 2012
# DESCR.:  Apply a median filter in spectral direction to the MUSE Datacube 
#          and subtract the filter from the cube
#          (do not alter the variances though, i.e. median filter is only applied
#          to flux datacube)

__version__ = "1.0.2"

import sys
import os
import argparse
import multiprocessing
import time
import pylab as p
from astropy.io import fits
from scipy.ndimage import filters
import line_em_funcs as lef

# get string of the commandline that was entered by the user
command = os.path.basename(sys.argv[0])
for entity in sys.argv[1:]:
    command = command+' '+entity

starttime = time.time()
get_timestring = lef.get_timestring

def medfilt_spectrum(data_part,filt_width):   
    """
    filtered_data_part = medfilt_spectrum(data_part,filt_width)
    the loop which does the actual median filtering on the alligned
    array of spectra
    """
    length = data_part.shape[1]
    medfilt_data = p.zeros_like(data_part)
    for i in xrange(length):
        # scipy.signal.medfilt is unbeliveable slow,
        # medfilt_data[:,i] = signal.medfilt(data_part[:,i],filt_width)
        medfilt_data[:,i] = filters.median_filter(data_part[:,i],
                                                  size=filt_width,
                                                  mode='mirror')
    return medfilt_data

# http://docs.python.org/library/argparse.html
parser = argparse.ArgumentParser(description="""Subtract an in spectral direction median-filtered version of the
 datacube.  Can be used to remove sources that have significant
 detectable continuum signal within the datacube.  """)
parser.add_argument("fitscube",
                    type=str,
                    help="The FITS File that contains the flux- and variances "
                         "datacube")
parser.add_argument("-S","--signalHDU",default=1,type=int,
                    help="HDU number (0-indexed) of fitscube, that contains the flux "
                    "(default: 1).")
parser.add_argument("-V","--varHDU",default=2,type=int,
                    help="HDU number (0-indexed) number of fitscube, that contains the variances."
                    "(default: 2).")
parser.add_argument("-o","--output",default=None,
                    help="Name of the output FITS file. "
                    "(default: median_filtered_W+<width>+<fitscube>)")
parser.add_argument("-W","--width",default=151,type=int,
                    help="The full width of the median filter. "
                    "If an even number is given, the actual filter value will "
                    "be this number +1. (Default: 151)")
parser.add_argument("-t","--num_cpu",default=None,type=int,
                    help="Number of CPUs which should be used "
                    "(default: all available).")
parser.add_argument("--memmap",action='store_true',
                    help="Use memory mapping. Reduces memory usage, but might also reduce speed of execution.)")
args = parser.parse_args()

# number of threads, either defined by user or use all CPUs
num_threads =  args.num_cpu
if num_threads == None:
    num_threads = multiprocessing.cpu_count()

width = args.width
if width % 2 == 0:
    # filter width must be uneven
    width = width + 1

fitscube = args.fitscube
signalHDU = args.signalHDU  # hdu indices are passed to the programm in 0-indexed
varHDU = args.varHDU        # notation
output = args.output

memmap = True if args.memmap else False

if output == None:
    # default output file name
    output = 'median_filtered_W'+str(width)+'_'+fitscube

print(fitscube+': Reading input flux (HDU'+str(signalHDU)+')... '+\
          get_timestring(starttime))
fitscube_hdu_list = fits.open(fitscube)
data, header = lef.read_hdu(fitscube,signalHDU, memmap=memmap)

# reshape the data for 1D-iteration (..so we only need no slow nested-loop)
print(str(fitscube)+': Prepare the flux data for spectral median filtering... '+\
          get_timestring(starttime))
shape_0 = data.shape[0]
shape_1 = data.shape[1]
shape_2 = data.shape[2]
num_spec = shape_1*shape_2
data = data.reshape(data.shape[0],shape_1*shape_2)

# multiprocessing median filtering starts
pool = multiprocessing.Pool(processes=num_threads)
threads = []
for j in range(num_threads):
    start = j*(num_spec/num_threads)
    end = (j+1)*(num_spec/num_threads)
    endstr = str(end)
    if j+1 == num_threads:
        end = num_spec + 1
        endstr = str(num_spec)
    print(str(fitscube)+': Core #'+str(j+1)+': Median Filtering spectra '+
          str(start+1)+'-'+endstr)
    data_part = data[:,start:end]
    threads.append(pool.apply_async(medfilt_spectrum,
                                    args=(data_part,width)))

pool.close() # no more workers are needed in the pool now

result = []
for t in threads:
    result.append(t.get()) # store the workers results in list result

pool.join() # wait until all threads are finished

print(str(fitscube)+': Writing out the median filtered datacube... '+\
          get_timestring(starttime))

# free some memory
# del data
del threads

medfiltered_data = p.concatenate(result,axis=1)
del result
# reshape to cube:
medfiltered_data = medfiltered_data.reshape(shape_0,
                                            shape_1,
                                            shape_2)
#print shape_0,shape_1,shape_2
# final data that is needed for object detection
med_sub_data = data.reshape(shape_0,shape_1,shape_2) - medfiltered_data

# prepare for writing output to disk
header['EXTNAME'] = 'MEDFILTERED_DATA'
header['MFR'] = (os.path.basename(sys.argv[0]),
                 'median filter subtract routine')
header['MFRC'] = (command,
                  'MFR full command called')
header['MFRV'] = (__version__,
                  'MFR version')
header['MFRW'] = (width,
                  'MFR median filter width')
header['MFRIN'] = (fitscube,
                   'ENR input FITS file')
header['MFRINS'] = (signalHDU,
                    'ENRIN flux HDU number')
header['MFRINN'] = (varHDU,
                    'ENRIN variance HDU number')
med_sub_header = header.copy()
med_sub_header['EXTNAME']  ='MEDFILTER-SUBTRACTED_DATA'
variance,varhead = lef.read_hdu(fitscube,varHDU,memmap=memmap)

# write final output for object detection
med_sub_data_hdu = fits.PrimaryHDU(data = med_sub_data,
                                   header = med_sub_header)
medfiltered_data_hdu = fits.ImageHDU(data = medfiltered_data,
                                     header = header)
variance_hdu = fits.ImageHDU(data = variance,
                             header = varhead)
# input cube contains exposure count HDU
#try:
#    ex_hdu = fitscube_hdu_list['EXP']
#    out_hdu_list = fits.HDUList(hdus=[med_sub_data_hdu,
                                      #variance_hdu,
                                      #medfiltered_data_hdu,
                                      #ex_hdu])
#except:
out_hdu_list = fits.HDUList(hdus=[med_sub_data_hdu,variance_hdu])

    
out_hdu_list.writeto(str(output),clobber=True,output_verify='silentfix')
print(str(fitscube)+
      ': All done! Filtered datacube (with unchanged variances) stored in: '+
      str(output)+
      ' (HDU 1: Original minus Medianf., HDU2: Variances, HDU3: Medianf.) - '+
      get_timestring(starttime))
