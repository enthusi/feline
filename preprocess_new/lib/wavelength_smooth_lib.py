# FILE:   wavelength_smooth_lib.py
# DESCR.: LSDcat Library for spectral cross-correlation
# AUTHOR: E.C.Herenz - eherenz@eso.org
# YEAR: 2020, 2023
# Copyright 2020, 2023 Edmund Christian Herenz
# LICENSE: BSD 3-Clause License - https://opensource.org/licenses/BSD-3-Clause

# FOR DETAILED DESCRIPTION OF THE STRUCTURE CREATED IN
# create_filter_matrix AND create_filter_matrix_vel
# see ../doc/spec_filt_notes/spec_filt_notes.pdf

# References:
# - Herenz, E. C. 2023, Astron. Nachr., e909.
# - Herenz, E. C. & Wisotzki, L. 2017, A&A, A111.

import math as m
import numpy as np
from scipy import signal
from scipy.sparse import csr_matrix
import getopt
import sys
import line_em_funcs as lef  # my own library with convenience functions
import multiprocessing
from astropy import constants

import warnings
warnings.filterwarnings("ignore",category=FutureWarning)

#  macOS since Sierra uses "spawn" 
try:
    multiprocessing.set_start_method('fork')
except RuntimeError:
    pass

def create_filter_matrix(delta_lambda, lambda_start=4800, lambda_em=1215.67,
                         cdelt=1.25, lMax=3463, trunc_constant=4, var_spec=None):
    """Creates the filter matrix (Eq. 15 in spec_filt_notes.pdf)

    Parameters
    ----------
    delta_lambda : float
        FWHM (in Angstrom / rest-frame) of Gaussian emission line
        template.
    lambda_start : float
        Wavelength (in Angstrom) of the first spectral pixel in the
        spectrum to be filtered.
    lambda_em : float
        Rest frame wavelength (in Angstrom) for the Gaussian line template.
        Default: 1215.67 Angstrom (Hydrogen Lyman Alpha).
    cdelt: float
        Wavelength increment per spectral pixel (in Angstrom;
        usally stored in CDELT1 in 1D spectra or CD3_3 in IFS datacubes).
    lmax: int
        Number of spectral pixels (e.g., from NAXISx keyword).
    trunc_constant:  float
        Constant that controls the truncation length for the filter (see Eq. 13
        in spec_filt_notes.pdf; defaults = 4).
    var_spec: np.ndarray
        Variance spectrum for formal matched filter template (default: None).
        If None, then only a normalised Gaussian will be used as the template.

    Returns
    -------
    filter_matrix : np.ndarray
        The lower-banded tri-diagonal matrix G (Eq. 15 in spec_filt_notes.pdf)
        that, when multiplied with a column-vector f (length lMax) that 
        represents the spectrum, will deliver the matched filter output for the
        template.  (Note: For the multiplication it is benificial to convert
        the matrix into the CSR sparse matrix representtion.)

    Notes
    -----
        See also https://doi.org/10.5281/zenodo.6471629 for the difference
        when using this function with and without a variance spectrum.

    """
    
    # convert all values to floats:
    delta_lambda = float(delta_lambda); lambda_start=float(lambda_start); 
    lambda_em = float(lambda_em); cdelt = float(cdelt); lMax = float(lMax)
    sqrt2pi_reci = 1./m.sqrt(2*m.pi) # normalization cnst.
    
    # calculation of fixed truncation length of the Gaussian filter:
    fwhm_factor = 2*m.sqrt(2*m.log(2))
    C = trunc_constant # <- constant for the truncation formula, literature says
          # should be chosen between 3 and 6
    
    # truncate all filters as the biggest window,
    # i.e. at the end of the wavelength range if v=const.
    tMax = abs(
        ((delta_lambda/lambda_em)*((lambda_start/cdelt)+lMax))/fwhm_factor)
    M = C*m.sqrt(tMax)+1 # ...numerically sufficient (I hope..)

    # .. since the sum runs from -M to +M we have to create vectors of the length 2M
    M = int(m.ceil(2*M))
    # CORRECT START VALUE (lambda_start) FOR EDGE EFFECTS
    lambda_start = lambda_start - cdelt*(M/2)
    if M % 2 == 0: # symmetric gaussian, zero indexed  
        M = M + 1  # so the maximum gets sampled, if we have an even
                   # number of elements (odd digit, zero indexed)

    # initialising the array that will hold in each line the filter profile
    # for each wavelength bin 
    h_l = np.zeros((int(lMax) + M - 1, M))                                     

    # filling the array with the filter profiles
    for i in range(int(lMax)+M-1):
        t = abs(((delta_lambda / lambda_em) * \
                 ((lambda_start / cdelt) + i))/fwhm_factor)
        # TODO: I sense a small inaccuracy here ^^
        # it should be ((lambda_start / cdelt) + i - M + 1) in the second line
        # does not matter in practice as t is linearly varying with i,
        # and M << lmax ... its also not correctly written in spec_filt_notes.pdf
        h_l[i,:] = (sqrt2pi_reci / t) * signal.windows.gaussian(M, t)
        
    # with h_l we can now fill the matrix for the
    # convolution as a dot product.
    filter_matrix = np.zeros((int(lMax) + M - 1, int(lMax)))
    for i in range(int(lMax)+M-1):
        # the matrix is a banded lower triangular matrix, so when
        # filling we have to check if we are in the upper left corner,
        # the central part or the lower right corner:
        if i < M: # upper left corner
            filter_matrix[i,0:i+1] = h_l[i,M-i-1:M] 
        elif i >= M and i < int(lMax): # central part
            filter_matrix[i,i+1-M:i+1] = h_l[i,:] 
        elif i >= int(lMax): # lower right corner
            filter_matrix[i,i+1-M:int(lMax)] = h_l[i,0:M-(i+1-int(lMax))]

    if var_spec is not None:
        assert len(var_sepec) == lMax, \
            "Variance spectrum has different length than expected."
        norm = 1. / np.sqrt(np.sum(filter_matrix**2 / var_spec, axis=1))
        filter_matrix /= var_spec
        filter_matrix = (filter_matrix.T * norm).T
        
    # NOTE: we have to deal with the edges, ...the main thing
    # is, that when multiplied to spectral vector, the resulting vector
    # will be longer -> so we truncate accordingly in the main routine
    # wavelength_smooth_parallel.py - but for now return the matrix as defined here

    return filter_matrix 


def create_filter_matrix_vel(velocity, lambda_start=4800., cdelt=1.3, lMax=3463,
                             trunc_constant=4., var_spec=None):
    """Creates the filter matrix for Gaussian line profiles; line width
       given as velocity (see spec_filt_notes.pdf).

    Parameters
    ----------
    velocity : float
        FWHM (in km/s - rest frame) of Gaussian emission line
        template.
    lambda_start : float
        Wavelength of the first spectral pixel in the spectrum that
        will be filtered.
    cdelt: float
        Wavelength increment per spectral pixel (in Angstrom;
        usally stored in CDELT1 in 1D spectra or CD3_3 in IFS datacubes).
    lMax: int
        Number of spectral pixels (e.g., from NAXISx keyword).
    trunc_constant:  float
        Constant that controls the truncation length for the filter (see Eq. 13
        in spec_filt_notes.pdf; defaults = 4).
    var_spec: np.ndarray
        Variance spectrum for formal matched filter template (default: None).
        If None, then only a normalised Gaussian will be used as the template.

    Returns
    -------
    filter_matrix : np.ndarray
        The lower-banded tri-diagonal matrix G (Eq. 15 in spec_filt_notes.pdf)
        that, when multiplied with a column-vector f (length lMax) that 
        represents the spectrum, will deliver the matched filter output for the
        template.  (Note: For the multiplication it is benificial to convert
        the matrix into the CSR sparse matrix representtion.)

    Notes
    -----
        See also https://doi.org/10.5281/zenodo.6471629 for the difference
        when using this function with and without a variance spectrum.

    """
    # same as create_filter_matrix but this time uses velocity as input - 
    # instead of \delta_\lambda and \lambda_0
    # to understand the source see create_filter_matrix function
    speedoflight =  constants.c.to('km/s').value 
    C = trunc_constant #truncation length
    sqrt2pi_reci = 1./m.sqrt(2*m.pi)
    fwhm_factor = 2*m.sqrt(2*m.log(2))
    velocity=float(velocity); lambda_start=float(lambda_start); 
    cdelt = float(cdelt); lMax = float(lMax)

    # this is the the first line which differs from  create_filter_matrix 
    #  lambda / delta_lambda -> velocity / speedoflight
    tMax = abs(
        ((velocity / speedoflight) * ((lambda_start / cdelt) + lMax))\
            /fwhm_factor
              )

    M = C*m.sqrt(tMax)+1  # should be big enough...
    M = int(m.ceil(2*M))
    lambda_start = lambda_start - cdelt*(M/2)
    if M % 2 == 0:
        M = M + 1 
    h_l = np.zeros((int(lMax)+M-1,M))
    for i in range(int(lMax)+M-1):  
        t = abs(
            ((velocity / speedoflight) * ((lambda_start / cdelt) + i))\
                /fwhm_factor
               )
        # see note in create_filter_matrix - the index i in above seems not super accurate,
        # and I think it should be i - M + 1 (shift because of edge effects)
        h_l[i,:] = (sqrt2pi_reci/t)*signal.windows.gaussian(M,t)

    filter_matrix = np.zeros((int(lMax)+M-1,int(lMax)))

    for i in range(int(lMax)+M-1):
        if i < M:
            filter_matrix[i,0:i+1] = h_l[i,M-i-1:M] 
        elif i >= M and i < int(lMax):
            filter_matrix[i,i+1-M:i+1] = h_l[i,:] 
        elif i >= int(lMax): 
            filter_matrix[i,i+1-M:int(lMax)] = h_l[i,0:M-(i+1-int(lMax))]

    if var_spec is not None:
        assert len(var_spec) == lMax, \
            "Variance spectrum has different length than expected."

        # LSDCat 2.0 - filter according to Eq. (25) in Sect. 2.3 of
        # Herenz (2023, AN paper)
        norm_specs = np.divide(np.square(filter_matrix), var_spec,
                               where=np.isfinite(var_spec))
        norm_denom = np.sqrt(np.sum(norm_specs, axis=1))
        norm = np.divide(1., norm_denom, where=(norm_denom != 0.))
        filter_matrix = np.divide(filter_matrix, var_spec, where=np.isfinite(var_spec))
        filter_matrix = (filter_matrix.T * norm).T
            
    return filter_matrix 


# the loop which does the actual filtering on the array of spectra
def filter_spectrum(filter_matrix, data):   
    """Dot-multiply filter matrix with each spectrum-vector in an array of spectrum vectors

    The spectral convolution operation in LSDCat is implemented as a
    (sparse) matrix multiplication; see ../doc/spec_filt_notes.pdf.
    Essentially h = G · f, where f is a vector storing the spectrum, G
    is the matrix for the convlution operation and h is the convolved
    spectrum.  

    The function filter_spectrum_parallel parallelises this function
    embarissingly by equipartitioning the spectrum vector array into
    similar sized chunks for each thread.

    Parameters
    ----------
    filter_matrix : np.ndarray
        The matrix G - Eq. 15 in spec_filt_notes.pdf
    data : np.ndarray
        The array storing the spectrum vectors f (Eq. 14 in
        spec_filt_notes.pdf): data[:,i] is the i-th spectrum.

    Returns
    -------
    filtered_data : np.ndarray
        The array storing the convolved vectors h = G · f (Eqs. 16 & 17 
        in spec_filt_notes.pdf).
        Note, in our implementation the filtered spetra need to be
        truncated left and right according to Eq. 18 in
        spec_filt_notes.pdf.

    """
    length = data.shape[1]
    filtered_data = np.zeros((filter_matrix.shape[0],length))
    filter_matrix_sparse = csr_matrix(filter_matrix)
    for i in range(length):
        # convolution implemented as matrix multiplication see
        # see ../doc/spec_filt_notes/spec_filt_notes.pdf
        filtered_data[:,i] = filter_matrix_sparse.dot(data[:,i])
        
    return filtered_data  


# the function that parallizes the filtering
def filter_parallel(filter_matrix,data,num_threads,string="spectra",
                    filename='',nans_select=None):
    """
     filtered_data = filter_parallel(filter_matrix,data,num_threads)
    the function that parallizes the filtering
    """
    assert len(data.shape) == 2 # i.e. data is reshaped
    shape = data.shape

    # determine the number of spectra which get processed
    if np.all(nans_select == None):
        length = shape[1]
        data_to_filter = data
    else:
        nans_select_num = np.sum(nans_select)
        length = shape[1] - nans_select_num

    pool = multiprocessing.Pool(processes=num_threads)
    threads = []
    nans_select_part_list = []

    for j in range(num_threads):
        start = j * int(length / num_threads)
        end = (j + 1) * int(length / num_threads)
        endstr = str(end)
        if j+1 == num_threads:
            # the last worker has to do some more work if
            # mod(spaxels,num_threads) != 0
            end=length+1 
            endstr=str(length)
        if np.all(nans_select != None):
            spectra_part = data[:,~nans_select][:,start:end]
        else:
            spectra_part = data[:,start:end]

        if filename is not None:
            print(str(filename)+': Thread '+str(j+1)+\
                  ' Working on '+string+' from #'+str(start+1)+\
                  ' to #'+endstr)
            
        threads.append(pool.apply_async(filter_spectrum,
                                        args=(filter_matrix,spectra_part)))
            
    pool.close() # no more workers are needed in the pool now

    result = []
    for t in threads:
        result.append(t.get()) # store the workers results in list result

    pool.join() # wait until all threads are finished
    data_dtype = data.dtype
    del data
    del threads

    # making 2D list of all spectra after "sub-sets" have been collected
    # from the workers:
    filtered_data = np.concatenate(result,axis=1)
    del result # free memory

    if np.all(nans_select != None):
        # ignored spectra in the filtering process are set to 0
        filtered_data_full = np.zeros((filter_matrix.shape[0],
                                      shape[1]), dtype=data_dtype)
        filtered_data_full[:,~nans_select] = filtered_data
        filtered_data = filtered_data_full

    return filtered_data

