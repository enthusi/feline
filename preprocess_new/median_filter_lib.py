import time
import multiprocessing
import numpy as np
from scipy.ndimage import filters

from line_em_funcs import get_timestring

#  macOS since Sierra uses "spawn"
try:
    multiprocessing.set_start_method('fork')
except RuntimeError:
    pass

def notch_interpol(spec, notch_start_px, notch_end_px, notch_bound):
    """
    mod_spec = notch_interpol(spec, filt_width, notch_start_px, notch_end_px, notch_bound)

    In:
    ---
    spec ... a spectrum, typically a datacube spaxel, as a 1D array 
    notch_start_px ... first spectral bin of the notch-filter blocked region
    notch_end_px ... last spectral bin of the notch-filter blocked region
    notch_bound ... size of the window before notch_start_px and after notch_end_point
                    to determine the base-points via  median over this spectral window.

    Out:
    ----
    mod_spec ... modified spectrum, filled with interpolated values over the notch-filter
                 blocked spectral bins.
    """

    # base points
    blue_base = np.median(spec[notch_start_px - notch_bound:notch_start_px + 1])
    red_base = np.median(spec[notch_end_px:notch_end_px + notch_bound + 1])
    # linear increment per pixel
    inc = (red_base - blue_base)/(notch_end_px - notch_start_px)
    # the interpolated values for each pixel in the gap
    int_inset = blue_base + np.indices((notch_end_px - notch_start_px,))[0] * inc
    # insert the interpolated values into the gap
    spec[notch_start_px:notch_end_px] = int_inset

    return spec


def medfilt_spectrum(data_part, filt_width, gwidth=0., ao=False, 
                     notch_start_px=0, notch_end_px=0, notch_bound=0):   
    """filtered_data_part = medfilt_spectrum(data_part, filt_width, gwidth
                                             ao=False, notch_start_px=0, 
                                             notch_end_px=0, notch_bound=0)

    The loop which does the actual median filtering and afterward Gaussian 
    filtering on a list of spectra.  With ao=True
    the notch-filter gap will be interpolated with notch_interpol.

    data_part[i,j] is an array with spectra, where the first index i 
    corresponds spectral pixels
    and the second index j indexes the spectra.

    """
    length = data_part.shape[1]
    medfilt_data = np.zeros_like(data_part) 

    if ao:
        for i in range(length):
            interp_data_part = notch_interpol(data_part[:,i],
                                              notch_start_px, notch_end_px, notch_bound)
            medfilt_data[:,i] = filters.median_filter(interp_data_part,
                                                         size=filt_width,
                                                         mode='mirror')
    else:
        for i in range(length):
            medfilt_data[:,i] = filters.median_filter(data_part[:,i],
                                                         size=filt_width,
                                                         mode='mirror')
    assert gwidth >= 0.
    if gwidth > 0:
        for i in range(length):
            medfilt_data[:,i] = filters.gaussian_filter1d(medfilt_data[:,i],
                                                          sigma=gwidth,
                                                          mode='nearest')
            
    return medfilt_data


def medfilt_cube_para(data, width, gwidth=0., 
                      num_threads=multiprocessing.cpu_count(), 
                      ao=False, notch_start_px=0, notch_end_px=0, 
                      notch_bound=0, verbose=True, fitscube='',
                      starttime=time.time()):
    if verbose:
        print(str(fitscube)+': Prepare the flux data for spectral median and Gauss filtering... '+\
              get_timestring(starttime))

    # reshape the data for 1D-iteration (..so we only need no slow nested-loop)
    shape_0 = data.shape[0]
    shape_1 = data.shape[1]
    shape_2 = data.shape[2]
    num_spec = shape_1*shape_2
    data = data.reshape(data.shape[0],shape_1*shape_2)

    if verbose:
        gwidth_info = ' / Post-Gauss filter g = '+str(gwidth) if gwidth != 0 else ''
        print(str(fitscube)+': Spectral filtering with a filter of W=' + \
              str(width)+ gwidth_info + '... '+get_timestring(starttime))
        

    # multiprocessing median filtering starts
    pool = multiprocessing.Pool(processes=num_threads)
    threads = []
    for j in range(num_threads):
        start = j*int(num_spec/num_threads)
        end = (j+1)*int(num_spec/num_threads)
        endstr = str(end)
        if j+1 == num_threads:
            end = num_spec + 1
            endstr = str(num_spec)
        print(str(fitscube)+': Core #'+str(j+1)+': Median Filtering spectra '+
              str(start+1)+'-'+endstr)
        data_part = data[:,start:end]
        if ao:
            threads.append(pool.apply_async(medfilt_spectrum,
                                            args=(data_part, width),
                                            kwds={'ao': True,
                                                  'notch_start_px': int(notch_start_px),
                                                  'notch_end_px': int(notch_end_px),
                                                  'notch_bound': notch_bound,
                                                  'gwidth': gwidth}))
        else:
            threads.append(pool.apply_async(medfilt_spectrum,
                                            args=(data_part, width),
                                            kwds={'gwidth': gwidth}))

    pool.close() # no more workers are needed in the pool now

    result = []
    for t in threads:
        result.append(t.get()) # store the workers results in list result

    pool.join() # wait until all threads are finished

    gf_str = ' (and Gaussian smoothing)' if gwidth != 0 else ''
    print(str(fitscube)+': Done with the median filtering' + gf_str + ' of the' + \
              ' datacube... ' + get_timestring(starttime))

    # free some memory
    del threads

    medfiltered_data = np.concatenate(result, axis=1)
    del result
    # reshape to cube:
    medfiltered_data = medfiltered_data.reshape(shape_0,
                                                shape_1,
                                                shape_2)
    
    return medfiltered_data
