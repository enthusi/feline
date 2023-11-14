# FILE:   spatial_smooth-lib.py
# AUTHOR: Edmund Christian Herenz
# DESCR.: Library of functions for spatial convolution of integral
#         field specttroscopy datacube layers.  Part of LSDCat.
#         
# LICENSE: BSD 3-Clause License // https://opensource.org/licenses/BSD-3-Clause
#
# If you make use of this code in your research please cite:
# - Herenz, E. C., & Wisotzki, L. 2017,  A&A 602, A111.
#   https://doi.org/10.1051/0004-6361/201629507
# - Herenz, E. 2023, AN, e606
#   https://doi.org/10.1002/asna.20220091

import warnings
import sys
import multiprocessing
import math as m
import pylab as p
import numpy as np
from scipy import signal

warnings.simplefilter("ignore", np.ComplexWarning) 

#  macOS since Sierra uses "spawn"
try:
    multiprocessing.set_start_method('fork')
except RuntimeError:
    pass

def spatial_filter_parts(kernel,data,method='normal'):
    """
    filtered_data = filter_parts(kernel,data)

    Filter each layer of 3D Datacube <data> with 2D Kernel from
    list of 2d kernels <kernel>.
    There must be the same number of elements in the list <kernel>
    as there are layers in the cube <data>

    """
    assert len(data.shape) == 3 and len(kernel) == data.shape[0]
    filtered_data = p.empty_like(data)
    length = data.shape[0]
    for i in range(length):
        if method == 'normal':
            filtered_data[i,:,:] = signal.convolve2d(data[i,:,:],kernel[i],
                                                     mode='same',
                                                     boundary='fill',
                                                     fillvalue=0)
        elif method == 'fft':
            filtered_data[i,:,:] = signal.fftconvolve(data[i,:,:],kernel[i],
                                                     mode='same')
        # TODO: try using astropy.convolution.convolve_fft
        # http://docs.astropy.org/en/stable/api/astropy.convolution.convolve_fft.html#astropy.convolution.convolve_fft
            
    return filtered_data    
filter_parts = spatial_filter_parts

def spatial_filter_parallel(kernel,cube,num_threads,selMask=False,mask=[],
                            filename='testdata',method='normal'):
    """
    filtered_data = filter_parallel(kernel,data,num_threads,selMask,mask=[])
    Parallizes the filter_parts function by splitting up the
    datacube in num_threads subcubes. Additionally masks
    can be specified (selMask -> use mask)
    """
    length = cube.shape[0]
    pool = multiprocessing.Pool(processes=num_threads)
    async_results = []

    # start all the workers....
    for j in range(num_threads):
        part_fac = int(length/num_threads)  # int(x) == int(m.floor(x))
        start = j*part_fac
        end = (j+1)*part_fac
        if j+1 == num_threads:
            end = length+1 # so also the last layer gets convolved!

        if filename is not None:
            print(str(filename)+': Thread '+str(j+1)+\
                  ': Working on wavelength layers from #'+str(start+1)+\
                  ' to #'+str(end))

        cube_part = cube[start:end,:,:]
        kernel_part = kernel[start:end]
        async_results.append(pool.apply_async(filter_parts,
                                        args=(kernel_part,
                                              cube_part,method)))

    pool.close()
    result = []
    for r in async_results:
        result.append(r.get())

    pool.join()
    del cube
    filtered = p.concatenate(result)
    if selMask == True and mask != []: # we dont apply the mask to the variance
        print(str(filename)+\
                  ': Applying mask to filtered data...')
        filtered *= mask
    del result  
    
    return filtered

# allow for shorthand
filter_parallel = spatial_filter_parallel

def prepArray_fixed(size):
    """x,y,x0,y0 = prepArray_fixed(size)
    
    Similar to numpy.meshrgrid, but for a 2D square array of size
    'size'.  If an even size is provided its made odd to have a well
    defined centre.  Returns x,y - coordinate arrays (see
    numpy.meshgrid) - and x0, y0 - coordinates of the centre.

    """
    if size % 2 == 0:
        size += 1
    x = p.arange(0, size, 1, dtype=float)
    y = x[:,p.newaxis]
    x0 = y0 = size // 2 # origin of coordinate system (int division)
    return x,y,x0,y0


def prepArray(trunc_constant, width_parameter):
    """
    x,y,x0,y0 = prepArray(trunc_constant,width_parameter)

    Preparation of the output array for 2D filtering windows.

    width_parameter ... width parameter of the filter window
                        (eg. sigma for Gaussian or R for Moffat)
    trunc_constant ... determines the size of the window, 
                       actial size will be trunc_constant*width_parameter + 1

    returns x,y,x0,y0 
    x,y -> axes , x0,y0 -> origin
    """
    R = float(width_parameter); C = float(trunc_constant)
    # C defines where filter windows is truncated (i.e. dimension of the array)
    # filter_window == 0 <=> r > trunc_constant*width_parameter + 1
    # width_parameter = sigma for Gaussian - fwhm for Moffat
    M = C*R + 1
    size = int(m.ceil(2*M))
    return prepArray_fixed(size)  # x,y,x0,y0

    
####
# FILTER WINDOWS
####

def makeMoffat(R,beta=3.5,trunc_constant=4, size=0):
    """
    M = makeMoffat(R, beta=3.5, trunc_constant=4, size=0)
    
    Makes a two-dimensional circular Moffat (Moffat, 1969 - A&A
    3:455-461) window for filtering (using FWHM & beta as input)

    R, beta ... parameters of the Moffat profile
    trunc_constant ... since only lim_(r->inf) M = 0, the Moffat
                       needs to be truncated at a certain radius
                       this is done at trunc_constant*fwhm + 1
    size ... set trunc_constant=0 and size>0 to define a fixed size of the output arrays
    """
    
    R = float(R)
    beta = float(beta)
    FWHM = 2*m.sqrt(2**(1/beta)-1)*R

    if trunc_constant > 0:
        # array is truncated at trunc_constant*FWHM + 1
        x,y,x0,y0 = prepArray(trunc_constant, FWHM)
    elif size >= 0:
        # array is truncated at a fixed size
        x,y,x0,y0 = prepArray_fixed(size)

    # Calcluate the Moffat profile; note that sampling effects are
    # neglected.
    norm = (beta - 1) / (m.pi * m.pow(R,2))
    M = (1 + ( (x-x0)**2+(y-y0)**2 )/m.pow(R,2) )**(-beta)
    M *= norm

    return M    

def makeMoffat_FWHM(FWHM, beta=3.5, trunc_constant=4, size=0):
    """
    G = makeMoffat_FWHM(FWHM,beta=3.5,runc_constant=4)

    see makeMoffat -> here FWHM is used as input parameter, insted of R

    FWHM ... FWHM of the Moffat
    beta ... Beta parameter of the Moffat profile
    trunc_constant & size ... see help(makeMoffat)
    """
    FWHM = float(FWHM)
    beta = float(beta)
    R = FWHM/(2*m.sqrt(2**(1/beta)-1))
    M = makeMoffat(R, beta=beta, trunc_constant=trunc_constant, size=size)
    return M


def makeGaussian_sigma(sigma,trunc_constant=4):
    """ 
    G = makeGaussian(sigma,trunc_constant=4)

    Makes a two-dimensional gaussian window for filtering

    sigma = std_deviation of the gaussian (in pixel units)
    trunc_constant = defines, where the filter window is set to 0
    (i.e. dimensions of the matrix)
                     i.e. G == 0 for all values > trunc_constant * sigma + 1
    """
    sigma = float(sigma)
    x,y,x0,y0 = prepArray(trunc_constant,sigma)

    # calculation of the gaussian:
    norm = 1./((2*m.pi)*m.pow(sigma,2))
    G =  np.exp( -(( x - x0 )**2 + ( y - y0 )**2) / ( 2*sigma**2 ))
    G *= norm

    return G

def makeGaussian_noncirc(sigma_x,sigma_y,theta_deg,trunc_constant=4):
    """
    sigma_x, sigma_y in pixel
    theta_deg in degree
    (currently not used in lsdcat)
    """
    sigma_max = np.max([sigma_x, sigma_y])
    x,y,x0,y0 = prepArray(trunc_constant, sigma_max)
    theta = m.radians(theta_deg)
    a = m.cos(theta)**2/(2*sigma_x**2) + m.sin(theta)**2/(2*sigma_y**2)
    b = m.sin(2*theta)/(2*sigma_x**2) - m.sin(theta)/(2*sigma_y**2)
    c = m.sin(theta)**2/(2*sigma_x**2) + m.cos(theta)**2/(2*sigma_y**2)
    G = np.exp(-a*(x-x0)**2 - b*(x-x0)*(y-y0) - c*(y-y0)**2)
    norm = np.sum(G)
    G /= norm
    return G


def makeGaussian_FWHM(fwhm,trunc_constant=4):
    """
    G = makeGaussian_FWHM(FWHM,trunc_constant=4)

    see help(makeGaussian_sigma)
    here the FWHM instead of sigma is given as the width parameter
    """
    sigma = fwhm / ( 2*m.sqrt(2*m.log(2)) )
    G = makeGaussian_sigma(sigma,
                           trunc_constant=trunc_constant)
    return G

def comp_sigma_poly(lbda,p):
    """ 
    sigma = comp_sigma_poly(lbda,p)

    Compute array of PSF widths (1sigma) as a function of wavelength.

    In:
    ---
    lbda ... 1D array of wavelengths 
    p ... Polynomial coefficents. The polynomial evaluates to the
          PSF FWHM over the wavelengths in lbda.

    Out:
    ----
    sigma ... Array containing the seeing width (1sigma) for each
              corresponding wavelength in lbda

    Note:
    -----
    Order of polynomial coefficiants are evaluated is:
    p[0]*x**(N-1) + p[1]*x**(N-2) + ... + p[N-2]*x + p[N-1].

    """
    
    fwhm = np.polyval(p,lbda)
    assert ~np.any(fwhm <= 0), 'FWHM <= 0... Check polynomial coefficients!'
    sigma = fwhm / (2*np.sqrt(2*np.log(2)))
    
    return sigma


def makeGaussians_sigmas(sigmas,trunc_constant):
    """
    gaussians = makeGaussians_sigmas(sigmas,trunc_constant)

    makes a list of 2D gaussians (as array) <gaussians> provided an
    array with sigmas
    """
    gaussians = []
    for sigma in sigmas:
        gaussians.append(makeGaussian_sigma(sigma,
                                            trunc_constant=trunc_constant))
    return gaussians

def makeGaussians_fwhms(fwhms,trunc_constant):
    """
    see help(makeGaussians_sigmas) -> here with FWHM as input parameter
    ... (not used at the moment - just for reference) ...
    """
    sigmas = []
    for fwhm in fwhms:
        sigma = float(fwhm)/(2*np.sqrt(2*np.log(2)))
        sigmas.append(sigma)
    gaussians = makeGaussians_sigmas(sigmas,trunc_constant)
    return gaussians

def makeMoffats_fwhms(fwhms, betas, trunc_constant, size):
    """
    moffats = makeMoffats(fwhms,betas,trunc_constant)

    Reaturn a list of Moffat profiles (each entry in the list is a
    numpy array).
    fwhms ... a list of fwhms
    betas ... a list of betas - or a constant value
    trunc_constant, size ... see help(makeMoffat)
    """
    assert len(fwhms) == len(betas)
    moffats = [makeMoffat_FWHM(fwhm,
                               beta=beta,
                               trunc_constant=trunc_constant,
                               size=size)
               for fwhm,beta in zip(fwhms, betas)]
    return moffats   
        
def makeGaussians_poly(xax, p, pix_scale=0.2, trunc_constant=4, classic=True):
    """
    gaussians = makeGaussians_poly(xax,p,
                                   pix_scale=0.2,trunc_constant=4):
    xax - 1D array containing the wavelengths of all slices in the MUSE cube
    p - array containing the coefficients of the polynomial (in arcsec)
        (note: the polynomial describes the FWHM - NOT the sigmas)
    pix_scale & trunc_constant - see help(makeGaussian)
    """
    sigmas = comp_sigma_poly(xax,p)
    sigmas /=  pix_scale
    gaussians = makeGaussians_sigmas(sigmas, trunc_constant)
    if not classic:
        # 10.5281/zenodo.6471630, Eq. 4, right hand side; but there ^2
        # missing in denominator (typo)
        gaussians = [gaussian / np.sqrt(np.sum(gaussian**2))
                     for gaussian in gaussians]
    
    return gaussians

def makeMoffats_poly(xax, p, b, pix_scale=0.2, trunc_constant=4, size=0, classic=True):
    """
    moffats = makeMoffats_poly(xax,p,betas,
                     pix_scale=0.2,trunc_constant=4)

    xax -- 1D array containing the wavelengths of all slices in the 
           MUSE cube
    p   -- array containing the coefficients of the polynomial 
           (in arcsec)
           (note: the polynomial describes the FWHM - NOT the sigmas)
    b - array containing the coefficients of the polynomial for beta
    pix_scale & trunc_constant & size - see help(makeMoffat)
    """
    fwhms = np.polyval(p, xax)
    assert ~np.any(fwhms <= 0), 'FWHM <= 0... Check polynomial coefficients!'
    fwhms /= pix_scale
    betas = np.polyval(b, xax)
    assert ~np.any(betas <= 0), 'beta <= 0... Check polynomial coefficients!'
    moffats = makeMoffats_fwhms(fwhms, betas, trunc_constant, size)
    if not classic:
        # see note in makeGaussians_poly
        moffats = [moffat / np.sqrt(np.sum(moffat**2))
                   for moffat in moffats]
        
    return moffats
        

# DUMP OF POTENTIALLY (HOPEFULLY) UNUSED STUFF BELOW

# # the next two functions are from the QSim code by R. Bacon
# # /work1/durham/qsim/source/dast/sciobj/trunk/lib/seeing.py
# # description of the formulas used can be found on the MUSE Wiki
# # https://musewiki.aip.de/fsf-model (original source still missing TODO)
# def comp_fwhm(lbda, seeing, l0=22, am=1.0):
#   """ return seeing FWHM values - for the a ground layer turbulence seeing modell
#   lbda: array of wavelengths in nm
#   seeing; dimm seeing FWHM at 500 nm and at zenith in arcsec
#   l0: turbulence outerscale length (m)
#   am: airmass
#   return the fwhm array
#   """
#   # all magic numbers are from the wiki... don't blame me!
#   fwhm0 = seeing*(am**0.6) * (500./lbda)**0.2
#   r0 = 2.01e-4*lbda/fwhm0
#   fwhm = fwhm0*np.sqrt(1.-2.183*(r0/l0)**0.356)
#   return fwhm

# def comp_sigma(lbda, seeing, l0=22, am=1.0):
#   """ return seeing sigma values (for gaussians)
#   lbda: array of wavelength in nm(!!)
#   seeing; dimm seeing FWHM(!!) at 500 nm and at zenith in arcsec
#   l0: turbulence outerscale length (m)
#   am: airmass
#   return the width parameter array
#   """
#   fwhm = comp_fwhm(lbda, seeing, l0, am)
#   sigma = fwhm / (2*m.sqrt(2*m.log(2)))
#   return sigma

# def makeGaussians_glt(xax,dimm_seeing,
#                   air_mass=1.0,scale_length=22.,pix_scale=0.2,trunc_constant=4):
#     """
#     gaussians = makeGaussians_glt(xax,dimm_seeing,
#                                air_mass=1.0,scale_length=22,
#                                pix_scale=0.2,trunc_constant=4)

#     xax - 1D array containing the wavelengths of all slices in the MUSE cube
#     dimm_seeing - dimm seeing FWHM in arcsec
#     air_mass - 1/cos(z), where z is the zenith distance of observation
#     scale_length - wavefront outer scale length [m], 22m typical for Paranal
#     (function is using a ground layer turbulence modell, taken from the 
#     MUSE Wiki: https://musewiki.aip.de/fsf-model ) 

#     returns list of Gaussians <gaussians> filter windows, with FWHMs calulated
#     with the ground layer turbulence modell
#     """
#     #  calculate sigmas in arcsecond (1D array)
#     sigmas = comp_sigma(xax,dimm_seeing,l0=scale_length,am=air_mass)
#     # transform sigmas to pixscale
#     sigmas /= pix_scale
#     gaussians = makeGaussians_sigmas(sigmas,trunc_constant)
#     return gaussians


# def makeMoffats_glt(xax,dimm_seeing,betas,
#                      air_mass=1.0,scale_length=22.,
#                      pix_scale=0.2,trunc_constant=4):
#     """
#     moffats = makeMoffats_glt(xax,dimm_seeing,betas,
#                      air_mass=1.0,scale_length=22.,
#                      pix_scale=0.2,trunc_constant=4)

#     same as makeGaussians_glt - for complete description see
#     help(makeGaussians_glt - here with moffat profiles).
#     betas ... either a float or a list with the beta value(s) 

#     returns list of Moffats <gaussians> filter windows, with FWHMs calulated
#     with the ground layer turbulence modell.
#     """
#     assert type(betas) == type(list()) or type(betas) == type(float())
    
#     # calcualte fwhms using the ground layer turbulence modell
#     fwhms = comp_fwhm(xax,dimm_seeing,
#                       l0=scale_length,am=air_mass)
#     # transform fwhms to pix scale
#     fwhms /= pix_scale
    
#     moffats = makeMoffats_fwhms(fwhms,betas,trunc_constant)
#     return moffats
    
