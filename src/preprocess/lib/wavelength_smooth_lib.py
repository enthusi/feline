# FILE:   wavelength_smooth_lib.py
# DESCR.: Library for functions used booth in wavelength_smooth.py and 
#         wavelength_smooth_parallel.py
# AUTHOR: C. Herenz (cherenz <at> aip.de)

# FOR DETAILED DESCRIPTION OF THE STRUCTURE CREATED IN
# create_filter_matrix AND create_filter_matrix_vel SEE MY
# MASTERTHESIS Appendix A.2
# http://knusper.net/master-thesis/Herenz-master.pdf

import math as m
import pylab as p
from scipy import signal
from scipy.sparse import csr_matrix
import getopt
import sys
import line_em_funcs as lef # my own library with convenience functions
import multiprocessing

import warnings
warnings.filterwarnings("ignore",category=FutureWarning)

def create_filter_matrix(delta_lambda,lambda_start=4800,lambda_em=1215.67,
                         cdelt=1.25,lMax=3463):

    """
    filter_matrix = create_filter_matrix(delta_lambda,
                     lambda_start=4800,lambda_em=1216,
                     cdelt=1.3,lMax=3463)

    ---
    Creates the filter matrix which will then be multiplied to
    a vector containing the spectrum.
    <delta_lambda> = linewidth of emission line (FWHM, in Angstrom)
                     for which the filter will be optimized
    <lambda_start> =  wavelength corresponding to 1st spectral pixel (in Angstrom)
                      (default: 4800 Angstrom)
    <lambda_em> = wavelength (in Angstrom) of emission line for which the filter 
                  will be optimized (default: Ly-alpha = 1216 Angstrom)
    <cdelt> = wavelength interval corresponding to each pixel (in Angstrom)
              (stored in CD3_3 in QSim Datacube-Header), default 1.3 Angstrom
    <lMax> = No. of spectral elements (stored in NAXIS3 in Qsim Datacube-Header)
    """
    # convert all values to floats:
    delta_lambda = float(delta_lambda); lambda_start=float(lambda_start); 
    lambda_em = float(lambda_em); cdelt = float(cdelt); lMax = float(lMax)

    # calculation of fixed truncation length of the Gaussian filter:
    fwhm_factor = 2*m.sqrt(2*m.log(2))
    C = 4 # <- constant for the truncation formula, literature says
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

    h_l = p.zeros((int(lMax)+M-1,M)) # empty matrix, will be filled
                                     # with all the different filters
                                     # (Gaussians at the moment - but
                                     # this can easily be changed)
                                     # every line of h_l will get its
                                     # Gaussian

    sqrt2pi_reci = 1./m.sqrt(2*m.pi) # normalization cnst.
    for i in xrange(int(lMax)+M-1):  
        t = abs(((delta_lambda/lambda_em)*((lambda_start/cdelt)+i))/fwhm_factor)
        h_l[i,:] = (sqrt2pi_reci/t)*signal.gaussian(M,t)

    # with h_l we can now fill the matrix for the
    # convolution as a dot product.
    filter_matrix = p.zeros((int(lMax)+M-1,int(lMax)))
    for i in xrange(int(lMax)+M-1):
        # the matrix is a banded lower triangular matrix, so when
        # filling we have to check if we are in the upper left corner,
        # the central part or the lower right corner:
        if i < M: # upper left corner
            filter_matrix[i,0:i+1] = h_l[i,M-i-1:M] 
        elif i >= M and i < int(lMax): # central part
            filter_matrix[i,i+1-M:i+1] = h_l[i,:] 
        elif i >= int(lMax): # lower right corner
            filter_matrix[i,i+1-M:int(lMax)] = h_l[i,0:M-(i+1-int(lMax))]

    # NOTE: we have to deal with the edges, ...the main thing
    # is, that when multiplied to spectral vector, the resulting vector
    # will be longer -> so we truncate accordingly in the main routine
    # wavelength_smooth_parallel.py - but for now return the matrix as defined here

    return filter_matrix 

def create_filter_matrix_vel(velocity,lambda_start=4800,cdelt=1.3,lMax=3463):
    """
    filter_matrix = create_filter_matrix_vel(velocity,
                                             lambda_start=4800,
                                             cdelt=1.3,lMax=3463)
    ---

    Creates the filter matrix which will then be multiplied to
    a vector containing the spectrum.

    <velocity> = linewidth in velocity space (i.e. constant rest-frame)
                 of emission line (FWHM, in km/s) for which the filter
                 will be optimized
    <lambda_start> =  wavelength corresponding to 1st spectral pixel (in Angstrom)
                      (default: 4800 Angstrom)
    <cdelt> = wavelength interval corresponding to each pixel (in Angstrom)
              (stored in CD3_3 in QSim Datacube-Header), default 1.3 Angstrom
    <lMax> = No. of spectral elements (stored in NAXIS3 in Qsim Datacube-Header)
    """
    # same as create_filter_matrix but this time uses velocity as input - 
    # instead of \delta_\lambda and \lambda_0
    # to understand the source see create_filter_matrix function
    speedoflight = 299792.458 # speed of light in km/s (with
                              # sufficient accuracy...)
    C = 4 #truncation length
    sqrt2pi_reci = 1./m.sqrt(2*m.pi)
    fwhm_factor = 2*m.sqrt(2*m.log(2))
    velocity=float(velocity); lambda_start=float(lambda_start); 
    cdelt = float(cdelt); lMax = float(lMax)

    # this is the the first line which differs from  create_filter_matrix 
    #  lambda / delta_lambda -> velocity / speedoflight
    tMax = abs(
        ((velocity / speedoflight) * ((lambda_start / cdelt)+lMax))\
            /fwhm_factor
              )

    M = C*m.sqrt(tMax)+1  # should be big enough...
    M = int(m.ceil(2*M))
    lambda_start = lambda_start - cdelt*(M/2)
    if M % 2 == 0:
        M = M + 1 
    h_l = p.zeros((int(lMax)+M-1,M))
    for i in xrange(int(lMax)+M-1):  
        t = abs(
            ((velocity / speedoflight) * ((lambda_start / cdelt) + i))\
                /fwhm_factor
               )
        h_l[i,:] = (sqrt2pi_reci/t)*signal.gaussian(M,t)
    filter_matrix = p.zeros((int(lMax)+M-1,int(lMax)))
    for i in xrange(int(lMax)+M-1):
        if i < M:
            filter_matrix[i,0:i+1] = h_l[i,M-i-1:M] 
        elif i >= M and i < int(lMax):
            filter_matrix[i,i+1-M:i+1] = h_l[i,:] 
        elif i >= int(lMax): 
            filter_matrix[i,i+1-M:int(lMax)] = h_l[i,0:M-(i+1-int(lMax))]
    return filter_matrix 


# the loop which does the actual filtering on the array of spectra
def filter_spectrum(filter_matrix,data):   
    """
    filtered_data = filter_spectrum(filter_matrix,data):
    the loop which does the actual filtering on the array of spectra
    """
    length = data.shape[1]
    filtered_data = p.zeros((filter_matrix.shape[0],length))
    filter_matrix_sparse = csr_matrix(filter_matrix)
    for i in xrange(length):
        # filtering implemented as matrix multiplication
        # see my Masters-Thesis for more description 
        # (see URL at the beginning of this file)
        # filtered_data[:,i] = p.dot(filter_matrix,data[:,i])
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
    if nans_select == None:
        length = shape[1]
        data_to_filter = data
    else:
        nans_select_num = p.sum(nans_select)
        length = shape[1] - nans_select_num

    pool = multiprocessing.Pool(processes=num_threads)
    threads = []
    nans_select_part_list = []

    for j in xrange(num_threads):
        start = j * (length / num_threads)
        end = (j + 1) * (length / num_threads)
        endstr = str(end)
        if j+1 == num_threads:
            # the last worker has to do some more work if
            # mod(spaxels,num_threads) != 0
            end=length+1 
            endstr=str(length)
        if nans_select != None:
            spectra_part = data[:,~nans_select][:,start:end]
        else:
            spectra_part = data[:,start:end]

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
    filtered_data = p.concatenate(result,axis=1)
    del result # free memory

    if nans_select != None:
        # ignored spectra in the filtering process are set to 0
        filtered_data_full = p.zeros((filter_matrix.shape[0],
                                      shape[1]), dtype=data_dtype)
        filtered_data_full[:,~nans_select] = filtered_data
        filtered_data = filtered_data_full

    return filtered_data


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
# # DEPRECATED STUFF BELOW, DELETE WHEN SURE IT IS NOT USED ANYWHERE ELSE

# def set_vars(parallel=True):
#     """
#     Set the variables for wavelength smoothing:
#     delta_lambda = set_vars(parallel=True)
#     """
#     p = parallel
#     # defaults
#     num_threads= 10
#     lambda_em = 1216
#     delta_lambda = 0.81 
#     data_hdu = 0 
#     stat_hdu = 1
#     outSel = False
#     try:
#         options,args = getopt.getopt(sys.argv[1:], "hi:o:F:l:S:N:t:",
#                                      ["help","input=","output=","FWHM=","lambda0=","SHDU=","NHDU=","threads="])
#     except getopt.error, msg:
#         usage(parallel=p)
#         sys.exit(2)
#     for option, value in options:
#         if option in  ("-h", "--help"):
#             usage(parallel=p)
#             sys.exit(0)
#         if option in ("-i","--input"):
#             inputfile = str(value)
#             selFile = True
#         if option in ("-o","--output"):
#             outputfile = str(value)
#             outSel = True
#         if option in ("-F","--FWHM"):
#             delta_lambda = float(value)
#         if option in ("-l","--lambda0"):
#             lambda_em = float(value)
#         if option in ("-S","--SHDU"):
#             data_hdu = int(value)-1
#         if option in ("-N","--NHDU"):
#             stat_hdu = int(value)-1
#         if option in ("-t","--threads"):
#             num_threads = int(value)

#     if outSel == False:
#         outputfile = 'wavelength_smooth_'+str(inputfile)

#     print(str(inputfile)+': Reading Signal...')
#     data,header = lef.read_hdu(inputfile,data_hdu)
#     print(str(inputfile)+': Reading Noise...')
#     stat,stat_header = lef.read_hdu(inputfile,stat_hdu)


#     if parallel == True:
#         return num_threads,inputfile,outputfile,data,header,stat,stat_header,delta_lambda,lambda_em
#     else:
#         return inputfile,outputfile,data,header,stat,stat_header,delta_lambda,lambda_em


# def usage(parallel=True):
#     # print out usage message    
#     if parallel == True:
#         sys.stdout.write(
# """
# wavelength_smooth_parallel.py - Wavelength smoothing of all spaxels in the datacube
#                                 with a gaussian kernel. Operation is performed on signal
#                                 and noise HDU. 
#                                 This version uses multi-threading.

# USAGE:

# wavelength_smooth_parallel.py [-h] -i <Input> [-o <outputfile> -F <delta_lambda> -l <lambda_em>
#                                -S <SignalHDU> -N <NoiseHDU> -t <NumThreads>]

# Possible options:
# -h [--help]: prints this message
# -i [--input]: specify input file (mandatory argument)
# -l [--lambda0]: rest-frame wavelength of emission line (in Angstrom) for which the
#                 optimal filter will be generated (default: 1216 - i.e. Ly-alpha)
# -F [--FWHM]: specify the FWHM of the gaussian kernel in Angstrom (default: 0.81) 
#              for which the optimal filter will be generated at rest-frame wavelength
#              <lambda_em>
# -S [--SHDU]: number of HDU containing the data (default: 1)
# -N [--NHDU]: number of HDU containing the noise (default: 2) 
#              (noise needs to be given as variance)
# -t [--threads]: number of threads (default: 10)

# """
# )
#     else:
#         sys.stdout.write(
# """
# wavelength_smooth.py - Wavelength smoothing of all spaxels in the datacube
#                        with a gaussian kernel. Operation is performed on signal
#                        and noise HDU. 

# USAGE:

# wavelength_smooth.py [-h] -i <Input> [-o <outputfile> -F <delta_lambda> -l <lambda_em>
#                                -S <SignalHDU> -N <NoiseHDU> -t <NumThreads>]

# Possible options:
# -h [--help]: prints this message
# -i [--input]: specify input file (mandatory argument)
# -l [--lambda0]: rest-frame wavelength of emission line (in Angstrom) for which the
#                 optimal filter will be generated (default: 1216 - i.e. Ly-alpha)
# -F [--FWHM]: specify the FWHM of the gaussian kernel in Angstrom (default: 0.81) 
#              for which the optimal filter will be generated at rest-frame wavelength
#              <lambda_em>
# -S [--SHDU]: number of HDU containing the data (default: 1)
# -N [--NHDU]: number of HDU containing the noise (default: 2) 
#              (noise needs to be given as variance)

# """
# )


# def create_filter_matrix(delta_lambda,lambda_start=4800,lambda_em=1216,
#                          cdelt=1.3,lMax=3463):

#     """
#     filter_matrix = create_filter_matrix(delta_lambda,lambda_start=4800,lambda_em=1216,cdelt=1.3,lMax=3463)

#     ---
#     Creates the filter matrix which will then be multiplied to
#     a vector containing the spectrum.
#     <delta_lambda> = linewidth of enmission line (FWHM, in Angstrom)
#                      for which the filter will be optimized
#     <lambda_start> =  wavelengt corresponding to 1st spectral pixel (in Angstrom)
#                       (default: 4800 Angstrom)
#     <lambda_em> = wavelength (in Angstrom) of emission line for which the filter 
#                   will be optimized (default: Ly-alpha = 1216 Angstrom)
#     <cdelt> = wavelength intervall corresponding to each pixel (in Angstrom)
#               (stored in CD3_3 in QSim Datacube-Header), default 1.3 Angstrom
#     <lMax> = No. of spectral elements (stored in NAXIS3 in Qsim Datacube-Header)
#     """
#     # convert all values to floats:
#     delta_lambda = float(delta_lambda); lambda_start=float(lambda_start); 
#     lambda_em = float(lambda_em); cdelt = float(cdelt); lMax = float(lMax)

#     # calculation of fixed truncation length of the gaussian filter:
#     fwhm_factor = 2*m.sqrt(2*m.log(2))
#     C = 4 # <- constant for the truncation formula, literature says
#           # should be choosen between 3 and 6
    
#     # truncate all filters as the biggest window,
#     # i.e. at the end of the wavelength range if v=const.
#     tMax = abs(((delta_lambda/lambda_em)*((lambda_start/cdelt)+lMax))/fwhm_factor)
#     M = C*m.sqrt(tMax)+1 # ...numerically sufficient (I hope..)

#     # .. since the sum runs from -M to +M we have to create vectors of the length 2M
#     M = int(m.ceil(2*M))
#     # CORRECT START VALUE (lambda_start) FOR EDGE EFFECTS
#     lambda_start = lambda_start - cdelt*(M/2)
#     if M % 2 == 0: # symmetric gaussian, zero indexed  
#         M = M + 1  # so the maximum gets sampled, if we have an even
#                    # number of elements (odd digit, zero indexed)

#     h_l = p.zeros((int(lMax)+M-1,M)) # empty matrix, will be filled
#                                      # with all the different filters
#                                      # (gaussians at the moment - but
#                                      # this can easily be changed)
#                                      # every line of h_l will get its
#                                      # gaussian

#     sqrt2pi_reci = 1./m.sqrt(2*m.pi) # normalization cnst.
#     for i in xrange(int(lMax)+M-1):  
#         t = abs(((delta_lambda/lambda_em)*((lambda_start/cdelt)+i))/fwhm_factor)
#         h_l[i,:] = (sqrt2pi_reci/t)*signal.gaussian(M,t)

#     # with h_l we can now fill the matrix for the
#     # convolution as a dot product.
#     filter_matrix = p.zeros((int(lMax)+M-1,int(lMax)))
#     for i in xrange(int(lMax)+M-1):
#         # the matrix is a banded lower triangular matrix, so when
#         # filling we have to check if we are in the upper left corner,
#         # the central part or the lower right corner:
#         if i < M: # upper left corner
#             filter_matrix[i,0:i+1] = h_l[i,M-i-1:M] 
#         elif i >= M and i < int(lMax): # central part
#             filter_matrix[i,i+1-M:i+1] = h_l[i,:] 
#         elif i >= int(lMax): # lower right corner
#             filter_matrix[i,i+1-M:int(lMax)] = h_l[i,0:M-(i+1-int(lMax))]

#     # NOTE: we have to deal with the edges, ...the main thing
#     # is, that when multiplied to spectral vector, the resulting vector
#     # will be longer -> so we truncate accordingly in the main routine
#     # wavelength_smooth_parallel.py - but for now return the matrix as defined here

#     return filter_matrix 
