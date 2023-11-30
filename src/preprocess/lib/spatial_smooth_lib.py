# FILE:   spatial_smooth-lib.py
# DESCR.: Functions for spatial_smooth_parallel.py
#         Spatial filtering with correct error propagation in datacube

import warnings
import sys
import multiprocessing
import math as m
import pylab as p
import numpy as np
from scipy import signal
from scipy import version

scipy_version = version.version
assert scipy_version[1] >= 9 or scipy_version[0] >= 1
warnings.simplefilter("ignore", np.ComplexWarning) 

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
    for i in xrange(length):
        if method == 'normal':
            filtered_data[i,:,:] = signal.convolve2d(data[i,:,:],kernel[i],
                                                     mode='same',
                                                     boundary='fill',
                                                     fillvalue=0)
        elif method == 'fft':
            filtered_data[i,:,:] = signal.fftconvolve(data[i,:,:],kernel[i],
                                                     mode='same')
            
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
        part_fac = length/num_threads
        start = j*part_fac
        end = (j+1)*part_fac
        if j+1 == num_threads:
            end = length+1 # so also the last layer gets convolved!
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
filter_parallel = spatial_filter_parallel

# NOTE:
# trunc_constant - defines where filter windows is truncated
# filter_window == 0 <=> r > trunc_constant*width_parameter + 1
# width_parameter = sigma for Gaussian - fwhm for Moffat

def prepArray(trunc_constant,width_parameter):
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
    C = float(trunc_constant)
    R = float(width_parameter)
    M = C*R + 1  # < constant for truncation, see page 53 of my master thesis
    # print M
    size = int(m.ceil(2*M))
    if size % 2 == 0: size += 1
    x = p.arange(0, size, 1, dtype=float)
    y = x[:,p.newaxis]
    x0 = y0 = size // 2 # origin of coordinate system
    return x,y,x0,y0

####
# FILTER WINDOWS
####

def makeMoffat(R,beta=3.5,trunc_constant=4):
    """
    M = makeMoffat(R,beta=3.5,trunc_constant=4)
    
    Makes a two-dimensional circular Moffat (Moffat, 1969 - A&A
    3:455-461) window for filtering (using FWHM & beta as input)

    R, beta ... parameters of the Moffat profile
    trunc_constant ... since only lim_(r->inf) M = 0, the Moffat
                       needs to be truncated at a certain radius
                       this is done at trunc_constant*fwhm + 1
    """
    
    R = float(R)
    beta = float(beta)
    FWHM = 2*m.sqrt(2**(1/beta)-1)*R
    
    # array is truncated at trunc_constant*FWHM + 1
    x,y,x0,y0 = prepArray(trunc_constant,FWHM)

    norm = (beta - 1) / (m.pi * m.pow(R,2))
    M = (1 + ( (x-x0)**2+(y-y0)**2 )/m.pow(R,2) )**(-beta)
    M *= norm

    return M    

def makeMoffat_FWHM(FWHM,beta=3.5,trunc_constant=4):
    """
    G = makeMoffat_FWHM(FWHM,beta=3.5,runc_constant=4)

    see makeMoffat -> here FWHM is used as input parameter, insted of R

    FWHM ... FWHM of the Moffat
    beta ... Beta parameter of the Moffat profile
    trunc_constant ... Moffat is trunccated at trunc_constant*fwhm + 1
    """
    FWHM = float(FWHM)
    beta = float(beta)
    R = FWHM/(2*m.sqrt(2**(1/beta)-1))
    M = makeMoffat(R,beta=beta,trunc_constant=trunc_constant)
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
    G =  p.exp( -(( x - x0 )**2 + ( y - y0 )**2) / ( 2*sigma**2 ))
    G *= norm

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

# the next two functions are from the QSim code by R. Bacon
# /work1/durham/qsim/source/dast/sciobj/trunk/lib/seeing.py
# description of the formulas used can be found on the MUSE Wiki
# https://musewiki.aip.de/fsf-model (original source still missing TODO)
def comp_fwhm(lbda, seeing, l0=22, am=1.0):
  """ return seeing FWHM values - for the a ground layer turbulence seeing modell
  lbda: array of wavelengths in nm
  seeing; dimm seeing FWHM at 500 nm and at zenith in arcsec
  l0: turbulence outerscale length (m)
  am: airmass
  return the fwhm array
  """
  # all magic numbers are from the wiki... don't blame me!
  fwhm0 = seeing*(am**0.6) * (500./lbda)**0.2
  r0 = 2.01e-4*lbda/fwhm0
  fwhm = fwhm0*np.sqrt(1.-2.183*(r0/l0)**0.356)
  return fwhm

def comp_sigma(lbda, seeing, l0=22, am=1.0):
  """ return seeing sigma values (for gaussians)
  lbda: array of wavelength in nm(!!)
  seeing; dimm seeing FWHM(!!) at 500 nm and at zenith in arcsec
  l0: turbulence outerscale length (m)
  am: airmass
  return the width parameter array
  """
  fwhm = comp_fwhm(lbda, seeing, l0, am)
  sigma = fwhm / (2*m.sqrt(2*m.log(2)))
  return sigma

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
    """
    
    fwhm = np.polyval(p,lbda)
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

def makeMoffats_fwhms(fwhms,betas,trunc_constant):
    """
    moffats = makeMoffats(fwhms,betas,trunc_constant)

    Reaturn a list of Moffat profiles (each entry in the list is a
    numpy array).
    fwhms ... a list of fwhms
    betas ... a list of betas - or a constant value
    trunc_constant ... truncation constant
    """
    moffats = []
    assert type(betas) == type(list()) or type(betas) == type(float())
    if type(betas) == type(list()):
        assert len(fwhms) == len(betas)
        for fwhm,beta in zip(fwhms,betas):
            moffats.append(makeMoffat_FWHM(fwhm,
                                           beta=beta,
                                           trunc_constant=trunc_constant))
    elif type(betas) == type(float()):
        beta = betas
        for fwhm in fwhms:
            moffats.append(makeMoffat_FWHM(fwhm,
                                           beta=beta,
                                           trunc_constant=trunc_constant))
    return moffats   
        
def makeGaussians_poly(xax,p,
                       pix_scale=0.2,trunc_constant=4):
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
    gaussians = makeGaussians_sigmas(sigmas,trunc_constant)
    return gaussians

def makeMoffats_poly(xax,p,betas,
                     pix_scale=0.2,trunc_constant=4):
    """
    moffats = makeMoffats_poly(xax,p,betas,
                     pix_scale=0.2,trunc_constant=4)

    xax -- 1D array containing the wavelengths of all slices in the 
           MUSE cube
    p   -- array containing the coefficients of the polynomial 
           (in arcsec)
           (note: the polynomial describes the FWHM - NOT the sigmas)
    betas - either a float or an list (same length as xax) with beta values
            for the Moffat profile.
    pix_scale & trunc_constant - see help(makeGaussian)                     
    """
    assert type(betas) == type(list()) or type(betas) == type(float())
    fwhms = np.polyval(p,xax)
    fwhms /= pix_scale
    moffats = makeMoffats_fwhms(fwhms,betas,trunc_constant)
    return moffats
        

def makeGaussians_glt(xax,dimm_seeing,
                  air_mass=1.0,scale_length=22.,pix_scale=0.2,trunc_constant=4):
    """
    gaussians = makeGaussians_glt(xax,dimm_seeing,
                               air_mass=1.0,scale_length=22,
                               pix_scale=0.2,trunc_constant=4)

    xax - 1D array containing the wavelengths of all slices in the MUSE cube
    dimm_seeing - dimm seeing FWHM in arcsec
    air_mass - 1/cos(z), where z is the zenith distance of observation
    scale_length - wavefront outer scale length [m], 22m typical for Paranal
    (function is using a ground layer turbulence modell, taken from the 
    MUSE Wiki: https://musewiki.aip.de/fsf-model ) 

    returns list of Gaussians <gaussians> filter windows, with FWHMs calulated
    with the ground layer turbulence modell
    """
    #  calculate sigmas in arcsecond (1D array)
    sigmas = comp_sigma(xax,dimm_seeing,l0=scale_length,am=air_mass)
    # transform sigmas to pixscale
    sigmas /= pix_scale
    gaussians = makeGaussians_sigmas(sigmas,trunc_constant)
    return gaussians

def makeMoffats_glt(xax,dimm_seeing,betas,
                     air_mass=1.0,scale_length=22.,
                     pix_scale=0.2,trunc_constant=4):
    """
    moffats = makeMoffats_glt(xax,dimm_seeing,betas,
                     air_mass=1.0,scale_length=22.,
                     pix_scale=0.2,trunc_constant=4)

    same as makeGaussians_glt - for complete description see
    help(makeGaussians_glt - here with moffat profiles).
    betas ... either a float or a list with the beta value(s) 

    returns list of Moffats <gaussians> filter windows, with FWHMs calulated
    with the ground layer turbulence modell.
    """
    assert type(betas) == type(list()) or type(betas) == type(float())
    
    # calcualte fwhms using the ground layer turbulence modell
    fwhms = comp_fwhm(xax,dimm_seeing,
                      l0=scale_length,am=air_mass)
    # transform fwhms to pix scale
    fwhms /= pix_scale
    
    moffats = makeMoffats_fwhms(fwhms,betas,trunc_constant)
    return moffats
    
