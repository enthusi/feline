# FILE: lsd_cat_lib.py
# DESCR: Routines for lsd_cat.py
# AUTHOR: C. Herenz (2013,2014,2015)

__version__ = "1.0.2"

import sys
import math as m
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from scipy.ndimage import measurements
#from astLib import astWCS   # http://astlib.sourceforge.net
                             # Replaced by astropy.wcs.WCS

import string

global numpy_version
global numpy_major_version
global numpy_minor_version

numpy_version = np.__version__.split('.')
numpy_major_version = int(numpy_version[0])
numpy_minor_version = int(numpy_version[1])

# global fits_version
# global fits_major_version
# global fits_minor_version

# fits_version = fits.__version__.split('.')
# fits_major_version = int(fits_version[0])
# fits_minor_version = int(fits_version[1])

# if fits_major_version == 3 and fits_minor_version >= 0:
#     import warnings
#     warnings.filterwarnings('ignore')

# check for numpy version - since in 1.6 a better counting
# algorithm for non_zero elements is present
assert numpy_major_version == 1
if numpy_major_version >= 6:
    def count_nonzero(array):
        return np.count_nonzero(array)
else:
    def count_nonzero(array):
        return np.sum(array != 0)

def tabvalue_test(tabvalues_list,out_line_form):
    for key in tabvalues_list:
        try:
            test = out_line_form[key]
        except KeyError:
            print('ERROR: TABVALUE KEY '+key+' NOT RECOGNIZED!'+\
                  ' ABORTING!')
            sys.exit(2)
    
    
def wavel(header,naxis=3,cdel_key='CD3_3'):
   """xax = wavel(header,naxis=3,cdel_key='CD3_3')

   header - a pyfits header object containing naxis, crval, crpix &
              crdel values needed for creation of wavelength grid
              array
   naxis - the wavelength axis (used to determine CRVALn, NAXISn &
             CRPIXn values)
   cdel_key - the key word that contains the wavelength increment
              array

   """

   nax = naxis

   naxis = header['NAXIS'+str(nax)]
   crval = header['CRVAL'+str(nax)]
   crpix = header['CRPIX'+str(nax)]
   crdel = header[cdel_key]

   xax = crval + (np.arange(naxis) - (crpix - 1))*crdel

   return xax
   


def pixelrad_from_cen(layer,x_cen,y_cen):
    """
    radii_to_cen = pixelrad_from_cen(image,x_cen,y_cen)

    Returns an array shaped like layer, with radial distances 
    from x_cen, y_cen in pixel coordinates.

    In:
    ---
    layer ... 2D array
    x_cen, y_cen ... coordinates of origin

    Out:
    ----
    radii_to_cen ... 2D array, with radial distances to origin

    """

    indices_y,indices_x = np.indices(layer.shape)  #  reverse axis-order!

    X = indices_x - x_cen
    Y = indices_y - y_cen

    radii_to_cen = np.sqrt(X**2+Y**2) 

    return radii_to_cen



def center_of_mass_weighted(weight_cube, label_cube, memory_friendly=False):
    """
    ids,x,y,z = center_of_mass_binary(labelcube,memory_friendly=False)

    Return center of mass of detection clusters (unweighted)
    ids,x,y,z,npix = center_of_mass_binary(labelczve)

    weight_cube
    labelcube = scipy.ndimage.measurements.label(sn_cube > thresh)[0]
    
    memory_friendly = True -> perform same calculation slower 
                              but more memory friendly (scales with no. of objs)
                              
    """

    # TODO: CURRENTLY THIS ROUTINE USES THE
    # SCIPY.MEASURMENTS.CENTER_OF_MASS ROUTINE - WHICH IS SLOW AND NOT VERY MEMORY EFFICIENT...
    # -> RECODE IN PURE NUMPY USING THE LABEL CUBE AND THE OBJECT SLICES

    assert label_cube.max() > 0 # TODO: does not work with no
                                # detections -- FIX!
    assert weight_cube.shape == label_cube.shape

    max_label = label_cube.max()
    ids = np.asarray(xrange(1,max_label+1))
    
    if memory_friendly == False:
        # using scipy.measurements center_of_mass routine for the calculation
        # - howver, it sucks quite some memory and is also not very fast
        center_of_masses = measurements.center_of_mass(weight_cube,
                                                       labels=label_cube,
                                                       index=xrange(1,max_label+1))
        # measurements returns a list of tuples (@X?!=wtf)
        center_of_masses = np.asarray(center_of_masses)[:,::-1]
        x_com = center_of_masses[:,0]
        y_com = center_of_masses[:,1]
        z_com = center_of_masses[:,2]

    else:
        # iterate over objects - slower, but more memory efficient
        x_com = np.empty(max_label)
        y_com = np.empty(max_label)
        z_com = np.empty(max_label)
        for label in xrange(1,max_label+1):
            center_of_mass = measurements.center_of_mass(weight_cube,
                                                         labels=label_cube,
                                                         index=label)

            x_com[label - 1] = center_of_mass[2]
            y_com[label - 1] = center_of_mass[1]
            z_com[label - 1] = center_of_mass[0]

    return ids,x_com,y_com,z_com


def center_of_mass_binary(label_cube,memory_friendly=False):
    """
    Return center of mass of detection clusters (unweighted)
    ids,x,y,z = center_of_mass_binary(labelcube,memory_friendly=False)

    labelcube = scipy.ndimage.measurements.label(sn_cube > thresh)[0]
    """
    com_un = center_of_mass_weighted(np.ones_like(label_cube),
                                     label_cube,
                                     memory_friendly = memory_friendly)
    return com_un


def calc_npix(label_cube,objects=None):
    """
    ids,npix = calc_npix(label_cube,objects=None)

    In:
    ---
    label_cube ... thresholded, labeled cube from 
                   scipy.measurements.label
    objects    ... a list of slices - one for the extent of each 
                   labeled object (output of 
                   scipy.measurments.find_objects(label_cube)
                   (optional - if not supplied, will be created 
                   on the fly, which might need some time)

    Out:
    ----
    ids        ... 1D numpy array with object ids 
    npix       ... 1D numpy array containging numbers of voxles 
                   belonging to detection
                   (i.e. ids[x] has npix[x] voxels)
    """
    max_label = label_cube.max()

    ids = xrange(1,max_label+1)
    if objects == None:
        objects = measurements.find_objects(label_cube) 

    npix = [count_nonzero(label_cube[obj_seg]) for obj_seg in objects]

    return np.asarray(ids),np.asarray(npix)


def calc_max(in_cube,label_cube,objects=None,windowed=False):
    """
    maxima = calc_max(in_cube,label_cube=None,objects=None)

    Calculate maximum value of in_cube for a detection.

    In:
    ---
    in_cube    ... cube - flux or detection significances (or something similar)
    label_cube ... cube with labeled detection clusters -
                   output of thresholding performed with
                   scipy.measurments.label
    objects    ... a list of slices - one for the extent of each 
                   labeled object (output of 
                   scipy.measurments.find_objects(label_cube)
                   (optional - if not supplied, will be created on the
                   fly, which needs some time)
    windowed   ... calculate only for detection cluster (i.e. similar to
                   windowed calculations on isophotoal contours in 2D)

    Out:
    ----
    maxima ... numpy.array of maximum values within detection clusters of in_cube

    """

    if objects == None:
        objects = measurements.find_objects(label_cube)

    if windowed == False:
        maxima = [np.max(in_cube[obj_seg]) for obj_seg in objects]
    else:
        maxima = [np.max(in_cube[obj_seg][obj_select]) for
                         obj_seg,obj_select in \
                             zip(objects, [np.nonzero(label_cube[obj_seg]) for
                                           obj_seg in objects])]
        

    return np.asarray(maxima)



def calc_min(in_cube, label_cube, objects=None, windowed=False):
    """
    minima = calc_det_sn_min(in_cube, label_cube, objects=None, windowed=False)
    
    same as calc_max -> but returns minima
    """
    if objects == None:
        objects = measurements.find_objects(label_cube)

    if windowed == False:
        minima = [np.min(in_cube[obj_seg]) for obj_seg in objects]
    else:
        minima = [np.min(in_cube[obj_seg][obj_select]) for
                         obj_seg,obj_select in \
                             zip(objects, [np.nonzero(label_cube[obj_seg]) for
                                           obj_seg in objects])]

    return np.asarray(minima)



def calc_median(in_cube, label_cube, objects=None, windowed=False):
    """
    medians = calc_median(in_cube, label_cube, objects=None, windowed=False)

    same as calc_max -> but returns median values in the detection clusters
    """
    if objects == None:
        objects = measurements.find_objects(label_cube)

    if windowed == False:
        medians = [np.median(in_cube[obj_seg]) for obj_seg in objects]
    else:
        medians = [np.median(in_cube[obj_seg][obj_select]) for
                   obj_seg,obj_select in \
                       zip(objects, [np.nonzero(label_cube[obj_seg]) for
                                     obj_seg in objects])]
        
        
    return np.asarray(medians)


def calc_stat(in_cube,label_cube,objects=None,windowed=False):
    """
    mean, stdev = calc_stat(in_cube, label_cube, objects=None) 

    Calculate mean & std.-deviation of in_cube for detection.

    -> see calc_max(..)

    Note: A tuple of 2 np.ndarrays is returned.
    """
    if objects == None:
        objects = measurements.find_objects(label_cube)

    if windowed == False:
        stat_seg = [[np.mean(in_cube[obj_seg]),
                     np.std(in_cube[obj_seg])] for obj_seg in objects]
    else:
        stat_seg = [[np.mean(in_cube[obj_seg][obj_select]),
                     np.std(in_cube[obj_seg][obj_select])] for
                    obj_seg,obj_select in \
                        zip(objects,
                            [np.nonzero(label_cube[obj_seg]) for
                             obj_seg in objects])]

    stat_seg = np.asarray(stat_seg)
 
    mean_seg = stat_seg[:,0]
    std_dev_seg = stat_seg[:,1]

    return mean_seg,std_dev_seg


def calc_npmax(in_cube, label_cube, peaks, pc=0.9, objects=None, windowed=False):
    """
    npmax = calc_npmax(in_cube, label_cube, peaks, pc=0.9,
                       objects=None, windowed=False)

    Calculate number of voxels that are above pc*peak within a detection.
    """
    assert pc > 0. and pc < 1.

    if objects == None:
        objects = measurements.find_objects(label_cube)

    assert len(objects) == len(peaks)

    if windowed == False:
        npmax = [count_nonzero(in_cube[obj_seg] > pc*peak) for
                 obj_seg, peak in zip(objects, peaks)]
    else:
        npmax = [count_nonzero(in_cube[obj_seg][obj_select] > pc*peak) for
                 obj_seg, obj_select, peak in zip(objects,
                                                  [np.nonzero(label_cube[obj_seg])
                                                   for obj_seg in objects],
                                                  peaks)]

    return np.asarray(npmax)


def calc_spikeyi(in_cube, label_cube, objects=None, windowed=False, r=0.5):
    """
    calc_spikeyi(in_cube, label_cube, objects=None, windowed=False, r=0.5):

    TODO: DOCUMENT 
    """

    if objects == None:
        objects = measurements.find_objects(label_cube)
    
    # sort (reverse) & flatten detections flux values
    if windowed == False:
        sorted_flattened = [np.sort(in_cube[obj_seg],axis=None)[::-1] 
                            for obj_seg in objects]
    else:
        sorted_flattened = [np.sort(in_cube[obj_seg][obj_select],
                                    axis=None)[::-1] 
                            for obj_seg,obj_select in
                            zip(objects,
                                [np.nonzero(label_cube[obj_seg]) for
                                obj_seg in objects])]

    # total sum
    summed = [np.sum(vals) for vals in sorted_flattened]

    # sorted cummulative sum
    cum_sum_sort = [np.cumsum(sorted_flat) for sorted_flat in sorted_flattened]

    # normed cum_sum
    cum_sum_normed = [cum_sum/summ for summ,cum_sum in zip(summed,
                                                          cum_sum_sort)]

    # find i (spikeyi parameter) where cum_sum_normed crosses i
    i_gt_r = [np.where(cum_sum_norm >= r)[0][0]
              for cum_sum_norm in cum_sum_normed]

    # upper and lower bounds around critical r
    r_up = [cum_sum_norm[i] if i > 1 else cum_sum_norm[1] 
            for i,cum_sum_norm in zip(i_gt_r,cum_sum_normed) ]
    r_down = [cum_sum_norm[i-1] if i-1 > 0 else cum_sum_norm[0] 
              for i,cum_sum_norm in zip(i_gt_r,cum_sum_normed)]

    # linear interpolation
    i_dash = [ (r-r_down_i)/(r_up_i-r_down_i) for r_up_i,r_down_i 
               in zip(r_up, r_down) ]

    # final spikyness ---> small = spikey ... larger = flat
    spikeyi = [i_dash_i + (i_gt_r_i - 1) for i_dash_i, i_gt_r_i in
               zip(i_dash,i_gt_r)]

    return np.asarray(spikeyi)
    

def calc_sum_pepi(signal_cube,label_cube,objects=None):
    """
    sums = calc_sum_pepi(signal_cube,label_cube)
    Sum values in minimal parallelepipeds in <signal_cube> definded around
    the object segments from <label_cube>. 
    """
    if objects == None:
        objects = measurements.find_objects(label_cube)
    sums = [np.sum(signal_cube[obj_seg]) for obj_seg in objects]
    return np.asarray(sums)


def calc_sum_seg(cube, label_cube, objects=None):
    """
    sum_seg = calc_sum_seg(det_sn_cube,label_cube,objects=None)
    calculate sum in the segments of the labels
    give a scipy.ndimage.measurements.find_objects output of label_cube 
    if at hand (otherwise it will be created)
    """
    if objects == None:
        objects = measurements.find_objects(label_cube)

    # nowlistening: DFTF (feat. DRS) by Need for Mirrors ;-)
    sum_seg = [np.sum(cube[obj_seg][obj_select]) for
               obj_seg,obj_select in \
                   zip(objects, [np.nonzero(label_cube[obj_seg]) for
                                 obj_seg in objects])]
    # :-) ^^ list-comprehension in a list-comprehension ^^ ;-)

    return np.asarray(sum_seg)


def calc_flux_in_aper(cube,x_cen,y_cen,z_cen,width,apradius,delta_lambda,
                      varcube=None):
    """flux(,vari)
    =
    calc_flux_in_aper(cube,x_cen,y_cen,z_cen,width,apradius,delta_lambda,
                      varcube=None)

    Summation of flux values in circular aperture over <width> layers
    centered on 'x_cen,y_cen,z_cen'.

    In:
    ---
    cube ... 3D fluxcube [erg/s/cm^2/A]
    x_cen,y_cen,z_cen ... position where aperture is centred
    width ... with of the window (in layers) of window over which the
              flux is summed
    apradius ... radius of the aperture (in pixels)
    delta_lambda ... wavelength increment (in Angstrom) per layer

    Out:
    ----
    flux ... flux within aperture [erg/s/cm^2]
    vari ... propagated variance within aper [erg/s/cm^2]**2 (only if varcube != None)

    """
    
    r_pix = pixelrad_from_cen(cube[0,:,:],x_cen,y_cen)
    app_sel = r_pix <= apradius

    narrow_band_cube = cube[z_cen - width/2.:z_cen + width/2.,
                            :,:]  
    narrow_band_image = np.sum(narrow_band_cube, axis=0)
    narrow_band_image *= delta_lambda  # erg/s/cm^2/A -> erg/s/cm^2

    flux = np.sum(narrow_band_image[app_sel])
    
    if varcube == None:
        return flux
    else:
        assert varcube.shape == cube.shape
        narrow_band_varcube = varcube[z_cen - width/2.:z_cen + width/2.,
                                      :,:]
        narrow_band_varimage = np.sum(narrow_band_varcube,axis=0)
        narrow_band_varimage *= delta_lambda**2
        variance = np.sum(narrow_band_varimage[app_sel])
        return flux,variance
    

# def calc_flux_in_aper_ellip(fluxcube,x_peak_sns,y_peak_sns,z_peak_sns,
#                             a_array,b_array,thetas_array,
#                             r_kron_array,delta_lambda_array,k=2.5,
#                             varcube=None):
#     """
#     flux(,vari)
#     =
#     calc_flux_in_aper_ellip(cube,x_cen,y_cen,z_cen,a,b,theta,r_kron,delta_lambda,
#                             varcube=None)

#     Summation of flux values in elliptical apperture defined via a, b, theta, r_kron.

#     In:
#     ---
#     fluxcube ... fluxcube
#     x_peak_sns,y_peak_sns,z_peak_sns ... peak SN coordinates of detections
#     a_array ... major axes of ellipses
#     b_array ... minor axes of ellipses
#     thetas ... position angle of ellipses (counterlockwise from north to east)
    

#     varcube (optional) ... variance cube

    
#     """

#     # thetas are in degrees, so convert to radians:
#     thetas = p.radians(thetas)

    

    

def calc_borders(label_cube,objects=None):
    """
    x_min,x_max,y_min,y_max,z_min,z_max 
    = 
    calc_borders(signal_cube,label_cube,objects=None)

    Calculates the corners of the minimum cube that contains a
    detection.
    """

    if objects == None:
        objects = measurements.find_objects(label_cube)

    coords = np.asarray([[obj_seg[0].start, obj_seg[0].stop,
                          obj_seg[1].start, obj_seg[1].stop,
                          obj_seg[2].start, obj_seg[2].stop] for obj_seg in objects])

    x_min = coords[:,4] ; y_min = coords[:,2] ; z_min = coords[:,0]
    x_max = coords[:,5] ; y_max = coords[:,3] ; z_max = coords[:,1]

    return x_min,x_max,y_min,y_max,z_min,z_max



def calc_maxima(cube, label_cube, objects=None):
    """
    maxima, max_x_coords, max_y_coords, max_z_coords
     =  
    calc_extrema(cube, label_cube,objects=None)

    Calculate extrema of detections in label_cube in cube.

    In:
    ---

    cube       ... cross-correlated S/N cube or flux cube
    label_cube ... thresholded, labeled cube from 
                   scipy.measurements.label
    objects    ... a list of slices - one for the extent of each 
                   labeled object (output of 
                   scipy.measurments.find_objects(label_cube)
                   (optional - if not supplied, will be created 
                   on the fly, which might need some time)

    Out:
    ----

    maxima       ... maximum value for each detection in cube
    max_x_coords ... maximum - x-coordinate (spatial)
    max_y_coords ... maximum - y-coordinate (spatial)
    max_z_coords ... maximum - z-coordinate (spectral)

    Note:
    -----

    This function is a replacement for what is provided by
    scipy.measurments, since the routines there were super slow and
    memory unefficent!

    """
    if objects == None:
        objects = measurements.find_objects(label_cube)

    maxima = calc_max(cube, label_cube, objects)
    minima = calc_min(cube, label_cube, objects)

    # borders of the subcubes bounding the detections
    x_min, x_max, y_min, y_max, z_min, z_max = \
                                    calc_borders(label_cube, objects)

    # maximum positions in the individual detection segments
    max_pos = np.asarray([np.unravel_index(cube[obj_seg].argmax(),
                                           cube[obj_seg].shape)
                          for obj_seg in objects])
    max_z_seg = max_pos[:,0]
    max_y_seg = max_pos[:,1]
    max_x_seg = max_pos[:,2]

    # transforming maximum positions to the coordinate system of the
    # full cube
    max_z_coords = z_min + max_z_seg 
    max_y_coords = y_min + max_y_seg
    max_x_coords = x_min + max_x_seg 

    return maxima, max_x_coords, max_y_coords, max_z_coords


def calc_borders_from_filter(x,y,z,cube,
                             wavel, wavel_unit,
                             ccs_p0,ccs_p1,ccs_p2,ccs_l0,ccl_vfwhm,
                             times_fwhm=3.,
                             spaxel_scale=0.2):
    """
    x,y,z - coordinates of detection 
    ccs_p0, ccs_p1, ccs_p2 - p0+p1*(l-l0)+p2*(l-l0)**2+p3*(l-l0)**3 coefficients
                             describing the wavelength dependence of 
                             PSF's FWHM
                             ([p0] - arcsec, [p1] - arcsec/AA,
                              [p2] - arcsec/AA^2 ...)
    ccs_l0 - wavelength zero point of the polynomial - i.e. p(l0) = p0
    times_fwhm - size of the window - centered on x,y,z
    ccl_vfwhm - v_FWHM of filter [km/s]
    spaxel_scale = 0.2 arcsec (default)
    """
    cube_dims = cube.shape
    cdelt = wavel[1] - wavel[0]
    
    x_dim = cube_dims[2]
    y_dim = cube_dims[1]
    z_dim = cube_dims[0]
    
    assert wavel_unit == 'Angstrom'  # later we might support
                                     # different units
    fwhm_arcsec = np.polyval([ccs_p2,ccs_p1,ccs_p0],
                             wavel - ccl_vfwhm)
    fwhm_xy_pix = fwhm_arcsec / spaxel_scale
    x_min = x - times_fwhm * fwhm_xy_pix[z.astype('int')]/2. 
    y_min = y - times_fwhm * fwhm_xy_pix[z.astype('int')]/2.
    x_max = x + times_fwhm * fwhm_xy_pix[z.astype('int')]/2.
    y_max = y + times_fwhm * fwhm_xy_pix[z.astype('int')]/2.

    fwhm_z_pix = ccl_vfwhm / 3e5 * wavel[z.astype('int')] / cdelt
   
    z_min = z - times_fwhm * fwhm_z_pix/2.
    z_max = z + times_fwhm * fwhm_z_pix/2.

    x_min[x_min < 0] = 0 
    y_min[y_min < 0] = 0 
    z_min[z_min < 0] = 0 

    x_max[x_max > x_dim] = x_dim 
    y_max[y_max > y_dim] = y_dim 
    z_max[z_max > z_dim] = z_dim 
    
    return x_min, x_max, y_min, y_max, z_min, z_max


def calc_weighted(x_min, x_max, y_min, y_max, z_min, z_max, cube,
                  thresh_ana=None, weigh_cube=None):
    """
    x_w, y_w, z_w = calc_weighted(x_min, x_max, y_min, y_max, z_min,
                                  z_max, cube, thresh_ana=None, weigh_cube=None)

    Descr.:
    -------
    Calculate windowed weighted coordinates.

    In:
    ---
    x_min, x_max, y_max, z_min, z_max ... arrays of boundary coordinates
                                          of mini-cubes (windows) in which
                                          calculation will be performed
    cube ... the sn_cube
    thresh_ana=None ... analysis threshold
    weigh_cube=None ... cube on which weighted coordinates above
                        above analysis threshold will be calculated (None means that cube is used.)

    Out:
    ----
    x_w, y_w, z_w ... arrays of weighted coordinates within the windows
                      in the cube coordinate system
    """

    x_w = [] ; y_w = [] ; z_w = []
    for x_min_i, y_min_i, z_min_i, x_max_i, y_max_i, z_max_i in \
        zip(x_min, y_min, z_min, x_max, y_max, z_max):
        # cut out subcube
        subcube = cube[int(z_min_i):int(z_max_i),
                       int(y_min_i):int(y_max_i),
                       int(x_min_i):int(x_max_i)]

        if weigh_cube is not None:
            weigh_subcube = weigh_cube[int(z_min_i):int(z_max_i),
                                       int(y_min_i):int(y_max_i),
                                       int(x_min_i):int(x_max_i)]

        if thresh_ana == None:
            if weigh_cube is None:
                com = measurements.center_of_mass(subcube)
            else:
                com = measurements.center_of_mass(weigh_subcube)
        else:
            assert thresh_ana >= 0
            if weigh_cube is None:
                subcube[subcube < thresh_ana] = 0
                com = measurements.center_of_mass(subcube)
            else:
                weigh_subcube[subcube < thresh_ana] = 0
                com = measurements.center_of_mass(weigh_subcube)
            
        x_com, y_com, z_com = np.asarray(com)[::-1]
        # translate back to original coordinate system
        x_w.append(x_com + x_min_i)
        y_w.append(y_com + y_min_i)
        z_w.append(z_com + z_min_i)

    return np.asarray(x_w), np.asarray(y_w), np.asarray(z_w)


def calc_max_win(x_min, x_max, y_min, y_max, z_min, z_max, cube):
    """
    x_peak, y_peak, z_peak = calc_max_win(x_min, x_max, y_min,
                                          y_max, z_min, z_max,
                                          cube)
    Descr.:
    -------
    Calculate windowed coordinate of maximum pixel in cube.

    In:
    ---
    x_min, x_max, y_max, z_min, z_max ... arrays of boundary coordinates
                                          of mini-cubes (windows) in which
                                          calculation will be performed
    cube ... the cube (e.g. fluxcube / sncube)

    Out:
    ----
    x_peak, y_peak, z_peak ... arrays of maximum coordinates within the windows
                               in the cube coordinate-system
    """

    x_peak = [] ; y_peak = [] ; z_peak = []
    for x_min_i, y_min_i, z_min_i, x_max_i, y_max_i, z_max_i in \
    zip(x_min, y_min, z_min, x_max, y_max, z_max):
        subcube = cube[z_min_i:z_max_i,
                       y_min_i:y_max_i,
                       x_min_i:x_max_i]

        peak_pos_i = np.unravel_index(subcube.argmax(), subcube.shape)
        
        z_peak.append(peak_pos_i[0] + z_min_i)
        y_peak.append(peak_pos_i[1] + y_min_i)
        x_peak.append(peak_pos_i[2] + x_min_i)
        
    return  np.asarray(x_peak), np.asarray(y_peak), np.asarray(z_peak)

        
def sn_2dmom(x_i,y_i,z_i,sn_cube,thresh_ana=3.0):
    """
    x_1mom,y_1mom,x_2mom,y_2mom,xy_2mom = \
       sn_2dmom(x_i,y_i,z_i,sn_cube,flux_cube,thresh_ana=3.0,
                 debug=False)

    Descr.:
    -------
    Calculate central moments in 2D for a detection. Preferably using 
    the layer of the SN-Peak of this detection (z_i).

    In:
    ---
    x_i,y_i,z_i ... SN peak coordinate of one detection
    sn_cube ... SN cube

    Out:
    ----
    x_1mom,y_1mom,x_2mom,y_2mom,xy_2mom 
       ... first and second central moments calculated in the layer
           where the SN peak

    """
        
    # 2D thresholding on peak layer with analysis threshold
    sn_peak_layer = sn_cube[z_i,:,:]
    sn_peak_logical = sn_peak_layer > thresh_ana
    sn_peak_label, sn_num_peaks = measurements.label(sn_peak_logical)

    # selecting the object corresponding to the inital detection
    obj_label = sn_peak_label[y_i,x_i]  # this is our guy...
    obj_isophotal_mask = sn_peak_label == obj_label
    obj_isophotal_vals = obj_isophotal_mask * sn_peak_layer

    # max_label=obj_label ensures that sn_peak_objects[-1] is the
    # slice for our guy
    sn_peak_objects = measurements.find_objects(sn_peak_label,
                                                max_label=obj_label)

    indices_y,indices_x = np.indices(obj_isophotal_vals.shape)

    # FIRST & SECOND IMAGE MOMENTS OF "ISOPHOTAL" SELECTED REGION
    x_1mom = np.sum(obj_isophotal_vals * indices_x) / \
             np.sum(obj_isophotal_vals)
    y_1mom = np.sum(obj_isophotal_vals * indices_y) /  \
             np.sum(obj_isophotal_vals)
    x_2mom = np.sum(obj_isophotal_vals * indices_x**2) / \
             np.sum(obj_isophotal_vals) - x_1mom**2
    y_2mom = np.sum(obj_isophotal_vals * indices_y**2) / \
             np.sum(obj_isophotal_vals) - y_1mom**2
    xy_2mom = np.sum(obj_isophotal_vals * indices_x * indices_y) / \
              np.sum(obj_isophotal_vals) - x_1mom * y_1mom


    return x_1mom,y_1mom,x_2mom,y_2mom,xy_2mom



def sn_2dmoms(x_peak_sns, y_peak_sns, z_peak_sns, sn_cube,
              thresh_ana=3.):
    """x_1moms, y_1moms, x_2moms, y_2moms, xy_2moms, sigma_isos
    sn_2dmoms(x_peak_sns, y_peak_sns, z_peak_sns, sn_cube, thresh_ana=3.)

    This is sn_2dmom for arrays of coordinates, and it also spills out
    sigma_isos.

    """
    x_1moms = [] ; y_1moms = [] ;  x_2moms = [] ;  y_2moms = []
    xy_2moms = [] ;  sigma_isos = []

    for x_peak_sns_i, y_peak_sns_i, z_peak_sns_i in zip(x_peak_sns,
                                                        y_peak_sns,
                                                        z_peak_sns):
        x_1mom_i, y_1mom_i, x_2mom_i, y_2mom_i, xy_2mom_i = \
            sn_2dmom(x_peak_sns_i, y_peak_sns_i, z_peak_sns_i,
                      sn_cube,thresh_ana=thresh_ana)
        sigma_iso_i = sigma_iso_circ(x_2mom_i, y_2mom_i)

        x_1moms.append(x_1mom_i)
        y_1moms.append(y_1mom_i)
        x_2moms.append(x_2mom_i)
        y_2moms.append(y_2mom_i)
        xy_2moms.append(xy_2mom_i)
        sigma_isos.append(sigma_iso_i)

    return np.asarray(x_1moms), np.asarray(y_1moms), np.asarray(x_2moms),\
        np.asarray(y_2moms), np.asarray(xy_2moms), np.asarray(sigma_isos)



def subcube_border_cuts(x_i,y_i,z_i, cube, ws=20):
    """
    xwin_min, xwin_max, ywin_min, ywin_max, zwin_min, zwin_max = \
         subcube_border_cuts(x_i,y_i,z_i, cube)
    """
    # cut out ~2ws x 2ws x 2ws subcube around x_i,y_i,z_i
    if x_i - ws <= 0:
        xwin_min = 0
    else:
        xwin_min = x_i - ws
    if x_i + ws >= cube.shape[2]:
        xwin_max = cube.shape[2]
    else:
        xwin_max = x_i + ws

    if y_i - ws <= 0:
        ywin_min = 0
    else:
        ywin_min = y_i - ws
    if y_i + ws >= cube.shape[1]:
        ywin_max = cube.shape[1]
    else:
        ywin_max = y_i + ws

    if z_i - ws <= 0:
        zwin_min = 0
    else:
        zwin_min = z_i - ws
    if z_i + ws >= cube.shape[0]:
        zwin_max = cube.shape[0]
    else:
        zwin_max = z_i + ws

    return int(xwin_min), int(xwin_max), int(ywin_min), int(ywin_max), int(zwin_min), int(zwin_max)



def cube_3dcom(x_i,y_i,z_i,sn_cube,thresh_ana=3.0,ws=20,
               anacube=None):
    """x_com, y_com, z_com = \
    cube_3dcom(x_i,y_i,z_i,sncube,thresh_ana=3.0,ws=20,
               anacube=None)

    3D center-of-mass (COM, in computer-vision paralance known as the
    first central moment) within an "isophotoal" area around the
    detection that is defined by the analysis threshold. Calculation
    is performed on a subcube that extends <ws> pixels in each
    direction from <x_i>,<y_i>,<z_i>.

    Note the difference to this libraries calc_weighted(...), which
    performs the calculation on a sub-cube (cuboid), without any any
    additional thresholding.

    Optionally an analysis cube can be supplied, in this case
    thresholding is performed on sncube, but the COM is calcualated on
    the analysis cube.

    In:
    ---
    x_i,y_i,z_i ... peak coordinate of detection
    sncube      ... SN Cube
    thresh_ana=3.0 ... analysis threshold
    ws=20. ... half maximum extension of 3D-subcube (in pix)
               to be analysed
    anacube ... the cube on which the COM will be computed

    Out:
    ----
    z_com,y_com,x_com ... COM coordinate (1st central moment in 3D)

    """

    xwin_min, xwin_max, ywin_min, ywin_max, zwin_min, zwin_max = \
                    subcube_border_cuts(x_i, y_i, z_i, sn_cube, ws=ws)
        
    subcube = sn_cube[zwin_min:zwin_max,
                     ywin_min:ywin_max,xwin_min:xwin_max]
    
    subcube_logical = subcube > thresh_ana
    subcube_logical_labels, subcube_num_det = \
                                measurements.label(subcube_logical)

    obj_label = subcube_logical_labels[z_i - zwin_min,
                                       y_i - ywin_min,
                                       x_i - xwin_min]

    if anacube == None:
        anasubcube = subcube
    else:
        anasubcube = anacube[zwin_min:zwin_max,
                             ywin_min:ywin_max,xwin_min:xwin_max]
    
    z_sub_com,y_sub_com,x_sub_com = \
            measurements.center_of_mass(anasubcube,
                                        labels=subcube_logical_labels,
                                        index=obj_label)

    z_com = z_sub_com + zwin_min
    y_com = y_sub_com + ywin_min
    x_com = x_sub_com + xwin_min

    return x_com, y_com, z_com


def cube_3dcoms(x_peak_sns, y_peak_sns, z_peak_sns, sn_cube,
                thresh_ana=3.0,anacube=None):
    """
    x_coms,y_coms,z_coms = \
           cube_3dcoms(x_peak_sns, y_peak_sns, z_peak_sns, sn_cube,
                       thresh_ana=3.0,anacube=None):

    see cube_3dcom for doc, here for arrays of peak coordinates
    """
    x_coms = [] ; y_coms = [] ; z_coms = []

    for x_i, y_i, z_i in zip(x_peak_sns, y_peak_sns, z_peak_sns):
        x_com_i, y_com_i, z_com_i = \
                cube_3dcom(x_i, y_i, z_i, sn_cube,
                           thresh_ana=thresh_ana,anacube=anacube)
        x_coms.append(x_com_i)
        y_coms.append(y_com_i)
        z_coms.append(z_com_i)

    return np.asarray(x_coms), np.asarray(y_coms), np.asarray(z_coms)


def sn_ana_flux(x_i,y_i,z_i, sn_cube, flux_cube, thresh_ana,
                varcube=None,ws=20):
    """
    flux(,vari) = sn_ana_flux(x_i,y_i,z_i, sn_cube, flux_cube, thresh_ana, ws=20)

    Summation of flux in analysis cluster above a certain SN threshold.

    In:
    ---
    x_i,y_i,z_i ... peak SN coordinate of a detection
    sn_cube ... S/N cube 
    flux_cube ... flux cube 
    var_cube ... variance cube
    thresh_ana ... analysis threshold to define analysis cluster via S/N cube
    varcube (optional) ... variance cube - if provided propagated variance will be
                           calculated
    ws ... +/- size of the window of the subcubes where analysis is performed
               (20 should be enough, i.e. a 40x40x40 cube centered on the emission line, 
               but for large objects (nearby galaxies or extended LyA) a bigger
               ws should be choosen)

    Out:
    ----
    flux ... sum of flux values in analysis cluster
    vari (only if variance cube is provided) ... propagated variance for flux 

    """
    xwin_min, xwin_max, ywin_min, ywin_max, zwin_min, zwin_max = \
                    subcube_border_cuts(x_i, y_i, z_i, sn_cube, ws=ws)

   
    flux_subcube = flux_cube[zwin_min:zwin_max,ywin_min:ywin_max,
                             xwin_min:xwin_max]
        
    sn_subcube = sn_cube[zwin_min:zwin_max,ywin_min:ywin_max,
                         xwin_min:xwin_max]

    if varcube != None:
        assert varcube.shape == flux_cube.shape
        varsubcube = varcube[zwin_min:zwin_max,ywin_min:ywin_max,
                             xwin_min:xwin_max]

    # get all voxels for object at x_i,y_i,z_i where sn > thresh_ana
    subcube_logical = sn_subcube > thresh_ana
    subcube_logical_labels, subcube_num_det = \
                                    measurements.label(subcube_logical)
    obj_label = subcube_logical_labels[z_i - zwin_min,
                                       y_i - ywin_min,
                                       x_i - xwin_min]
    select_cube = subcube_logical_labels == obj_label

    sel_im = p.sum(select_cube,axis=0)

    flux_vals = flux_subcube[select_cube]
    flux = p.sum(flux_vals)

    if varcube != None:
        vari_vals = varsubcube[select]
        vari = p.sum(vari_vals)
        return flux, vari

    else:
        return flux

    
    
def sn_ana_fluxes(x_peak_sns, y_peak_sns, z_peak_sns, sn_cube, flux_cube,
                  varcube=None,thresh_ana=3.0,ws=20.):
    """fluxes(,variances) 
    =
    sn_ana_fluxes(x_peak_sns, y_peak_sns, z_peak_sns, sn_cube,
    flux_cube, varcube=None,thresh_ana=3.0,ws=20.)

    
    Same as sn_ana_flux - but for arrays of (x,y,z)_peak_sns.
    """
    fluxes = []
    if varcube != None:
        variances = []
    
    for x_i, y_i, z_i in zip(x_peak_sns, y_peak_sns, z_peak_sns):
        if varcube != None:
            flux_i, vari_i = sn_ana_flux(x_i,y_i,z_i,
                                         sn_cube, flux_cube, thresh_ana,
                                         varcube=varcube, ws=ws)
            fluxes.append(flux_i)
            variances.append(vari_i)
        else:
            flux_i = sn_ana_flux(x_i,y_i,z_i, sn_cube, flux_cube, thresh_ana,
                                 varcube=None, ws=ws)
            fluxes.append(flux_i)

    if varcube != None:
        return p.asarray(fluxes), p.asarray(variances)
    else:
        return p.asarray(fluxes)




def sigma_iso_circ(x_2mom,y_2mom):
    """ sigma_iso = sigma_iso_circ(x_2mom,y_2mom):

    In:
    ---
    x_2mom, y_2mom - 2nd central moments in x and y direction

    Out:
    ----
    sigma_iso ... sqrt of 2nd central moment assuming circular 
                  symmetry
                  (std. deviation for a 2D circular Gaussian
                  distribution)

    """

    vari_iso = (x_2mom + y_2mom)/2.
    sigma_iso = np.sqrt( vari_iso )

    return sigma_iso


def ellipse_parameters(x_2mom,y_2mom,xy_2mom):
    """a,b,theta,elong,ellip  = ellipse_parameters(x2_mom,y_2mom,xy_2mom):

    Calculation of ellipse parameters from 2D 2nd central moments: 

    a ... major axis, b .... minor axis, elong ... elongation, ellip
    ... ellipticity

    (Following equations in Sect. 10.1.5 in SExtractor manual)

    In:
    ---
    x_2mom, y_2mom, xy_2mom - 2nd central moments (2D)

    Out:
    ----
    a,b,theta,elong,ellip ... ellipse describing the 2nd central moments
                              (theta in degrees)

    """
    
    vari_iso = (x_2mom + y_2mom)/2.

    aa = vari_iso + np.sqrt( vari_iso + xy_2mom**2 )
    bb = vari_iso - np.sqrt( vari_iso + xy_2mom**2 )
    a = np.sqrt(aa)
    b = np.sqrt(bb)

    theta = np.arctan(2 * (xy_2mom / (x_2mom - y_2mom)))

    # [-pi/2,pi/2[ ambiguty of arctan requires to 
    # theta needs to have the same sign as xy_2mom!
    # (arctan(-x) = - arctan(x))
    theta[np.sign(theta) != np.sign(xy_2mom)] *= -1  # I HOPE THIS IS CORRECT
    
    theta_deg = np.degrees(theta)

    elongation = a / b
    ellipticity = 1 - 1./elongation

    return a, b, theta_deg, elongation, ellipticity

    


def calc_lambda(z_coords,header):
    """
    lambdas = calc_lambda(z_coords,header)
    
    Calculates wavelengths for z-coordinates of datacube.

    header ... header of a datacube
    z_coords ... array containing (sub-)pixel coordinates  (zero indexed)
    returns: lambda(z_coords)
    """
    # assert type(header) == type(fits.header.Header())
    assert header['NAXIS'] == 3
    crpix = header['CRPIX3']; crval = header['CRVAL3']
    try:
        crdelt = header['CD3_3']
    except KeyError:
        try:
            crdelt = header['CDELT3']
        except KeyError:
            print('ERROR - INCOMPATIBLE FITS HEADER FOR WAVELENGTH CALCULATION!')
            sys.exit(2)
            
    lambdas = (z_coords-(crpix-1))*crdelt + crval

    return lambdas


def pix_to_radec(x_coords,y_coords,header): 
    """
    ra_coords, dec_coords = pix_to_radec(x_coords,y_coords,header)

    Calcualtes RA's & DEC's for x- & y-coordinates of datacube.

    In:
    ---
    header ... header of datacube, containing WCS information
    x_coords ... array containing the (sub-)x-pixelcoordinates (zero indexed)
    y_coords ... array containing the (sub-)y-pixelcoordinates (zero indexed)

    Out:
    ----
    ra,dec ... RAs & DECs corresponding to x_coords & y_coords
    
    """
	
    try :               #Possible to have input of a single x and y
                        #value pair
        size=len(x_coords)
    except TypeError:
        size=1

    #Create wcs obj
    wcs_obj=WCS(header)

    #Convert to (ra,dec)
    try:
        # this works only if header is from a 2D image
        (radec) = wcs_obj.wcs_pix2world(x_coords,y_coords,0)
    except TypeError:
        # this is needed if header is from a cube
        (radec) = wcs_obj.wcs_pix2world(x_coords,y_coords,np.zeros(size),0)

    if size !=1:
        return(radec[0],radec[1]) #radec is an array of 2 arrays : an
                                  #array containg all RAs, and an array
                                  #containing all DECs
    else:
        return([radec[0][0]],[radec[1][0]]) #If only one x,y pair, we
                                        #still write to a list to
                                        #prevent type errors


def gen_tempfilename(suffix='.fits',chars=6):
    """
    Generate temporary filename 'XXXXXXX.suffix', where 
    XXXXXX are <chars> (default: 6) random charchaters.
    """
    assert type(suffix) == str
    assert type(chars) == int

    tempfilename = ''.join(np.random.choice(list(string.ascii_uppercase + \
                                                 string.digits),size=chars))
    tempfilename += suffix
    return tempfilename


# CURRENTLY UNUSED FUNCTIONS - CHARACTERISTIC RADIUS
def calc_char_rad(signal_cube,
                  x_cen,y_cen,z_cen,
                  x_min,x_max,y_min,y_max,z_min,z_max,
                  char_method='char',
                  max_mult=1.):
    """
    char_rads,char_rads_flag =  calc_char_rad(signal_cube,
                                x_min,x_max,y_min,y_max,z_min,z_max):

    Calculates characteristic radii or compactness
    - in the style of Kron 1980
    - see Equation (3) or (4) in Infante, L. 1987, A&A, 183, 177 for
      definition of characteristic radius or compactness

    In:
    --
    signal_cube - the cube which is used to assign the weights for the
                  radius (should be either flux, or cross-correlated s/n)

    x_cen,y_cen,z_cen - central coordinates of the detections
    x_min,x_min,x_max,y_min,y_max,z_min,z_max - corner coordinates (arrays)
         of the minimum parallelepipeds that contain a detection 

    max_mult - by default the radius of the sphere in which the
               caracteristic radius will be determined is the maximum extend
               of the parallelepid - this can be altered by multiplication
               with this factor (default: 1. - i.e. original length)
    char_method - quantity to be calculated: 
                  'char' => characteristic radius
                  'com'  => compactness

    Out:
    ---
    char_rads - characteristic radii / or compactness 
                of the detections in signal cube
    char_rads_flag - the number of edges affecting the integration volume 
                     + the number of strange things happening during the calculation
                    (in general: the higher, the less trustworthy is the output -
                    its safe to trust those char_rads with char_rads_flag == 0)
    """
    assert char_method == 'char' or char_method == 'com'
    assert np.all(len(x_min) == len(x_max) == len(x_cen) == \
                  len(y_min) == len(y_max) == len(y_cen) == \
                  len(z_min) == len(z_max) == len(z_cen))

    num_dets = len(x_min)  # number of detections
    char_rads_flag = np.zeros(num_dets,dtype=int)  # flag values 
    z_cube_max,y_cube_max,x_cube_max = signal_cube.shape

    coordinate_tripel_list = [[x_min,x_max,x_cen],
                              [y_min,y_max,y_cen],
                              [z_min,z_max,z_cen]]

    # test if supplied coordinate tripels make sense - if not, flag them.
    for coordinate_tripel in coordinate_tripel_list:
        if not np.all((coordinate_tripel[0] <= coordinate_tripel[2]) \
                          & \
                      (coordinate_tripel[1] >= coordinate_tripel[0])):

            err_select_1 = coordinate_tripel[0] > coordinate_tripel[2]
            err_select_2 = coordinate_tripel[1] < coordinate_tripel[0]
            err_select = np.logical_or(err_select_1,err_select_2)

            char_rads_flag[err_select] += 1
            
    # upper limit of integration (r_up in Eq.(3) of Infante1987 - here
    # r_max)
    r_max = np.sqrt(  (x_max - x_min)**2 \
                    + (y_max - y_min)**2 \
                    + (z_max - z_min)**2 )

    r_max *= max_mult / 2. 
  
    # Calculations of r_char will be performd on subcubes.  These are
    # centered on each objects (x_cen,y_cen,z_cen). Boundaries of the
    # subcubes are tangent planes, that are parallel to the coordinate
    # axes and enclose the sphere of r_max.  Calculating coordinate
    # tripels (min + max) for cutting out those subcubes:
    subcube_boundary_list = [[np.floor(cen - r_max) - 1,
                              np.ceil(cen + r_max) + 1]
                              for cen in [z_cen,y_cen,x_cen]]
    [[z_subcube_min,z_subcube_max],
     [y_subcube_min,y_subcube_max],
     [x_subcube_min,x_subcube_max]] = subcube_boundary_list

    # For a detection near the datacube edges, parts of the r_max
    # sphere might lie outside of the cube. We ignore those voxels,
    # and flag those objects accordingly.
    for axes_min in [z_subcube_min,y_subcube_min,x_subcube_min]:
        min_sel = axes_min <= 0
        axes_min[min_sel] = 0
        char_rads_flag += min_sel.astype(np.int8)
    for ax_max_coord,axes_max in zip([z_cube_max,y_cube_max,x_cube_max],
                                     [z_subcube_max,
                                      y_subcube_max,
                                      x_subcube_max]):
        max_sel = axes_max >= ax_max_coord
        axes_max[max_sel] = ax_max_coord
        char_rads_flag += max_sel.astype(np.int8)
        
    # Center coordinates need also to be changed with respect to the
    # new subcube boundaries.
    subcube_center_list = [z_cen - z_subcube_min,
                           y_cen - y_subcube_min,
                           x_cen - x_subcube_min]
    z_subcube_cen = subcube_center_list[0]
    y_subcube_cen = subcube_center_list[1]
    x_subcube_cen = subcube_center_list[2]

    # DEBUG: check wether coordinates are really within sub-cubes
    for subcube_center_coordinate,max_coordinate in \
            zip(subcube_center_list,
                [z_subcube_max,y_subcube_max,x_subcube_max]):
        assert np.all(subcube_center_coordinate <= max_coordinate)

    
    # calculate the characterstic radius by integrating over r_max
    # sphere for every subcube:
    char_rads = np.zeros(num_dets,dtype=float)
    for i in xrange(num_dets):
        subcube = signal_cube[z_subcube_min[i]:z_subcube_max[i],
                              y_subcube_min[i]:y_subcube_max[i],
                              x_subcube_min[i]:x_subcube_max[i]]
        Z,Y,X = np.indices(subcube.shape,dtype=float)
        r_max_sphere =  np.sqrt((Z-z_subcube_cen[i])**2 + \
                                (Y-y_subcube_cen[i])**2 + \
                                (X-x_subcube_cen[i])**2)

        # define integration boundary by setting all
        # other values of subcube to 0
        r_max_sphere_select = r_max_sphere <= r_max[i]
        subcube_cut = subcube * r_max_sphere_select

        # actual integration
        if char_method == 'char':
            # Eq. 3 - Infante 1987
            r_eff_nom = np.sum(subcube_cut * r_max_sphere)  # nominator
            r_eff_denom = np.sum(subcube_cut)  # denominator
            r_eff = r_eff_nom / r_eff_denom  
        elif char_method == 'com':
            # Eq. 4 - Infante 1987
            r_eff_denom = np.sum(subcube_cut * (r_max_sphere**2)**(-1))
            r_eff_nom = np.sum(subcube_cut)
            if r_eff_nom < 0:  # integration might produce negative values - we discard those
                r_eff = -999.
                char_rads_flag[i] += 1
            elif r_eff_denom < 0:
                r_eff = -999.
                char_rads_flag[i] += 1
            else:
                r_eff = m.sqrt(r_eff_nom / r_eff_denom)
            
        # now - to be sure - also throw away NaNs
        if r_eff != r_eff:
            # only nans do not equal themselves
            r_eff = -999.
            char_rads_flag[i] += 1

        char_rads[i] = r_eff
        
    assert len(char_rads) == len(char_rads_flag)
    return char_rads, char_rads_flag


# def calc_sn_stat(det_sn_cube,label_cube,memory_friendly=False):
#     """
#     det_sn_mean,det_sn_std_dev = calc_sn_stat(det_sn_cube,
#                                               label_cube,
#                                               memory_friendly=False)
#     memory_friendly = True -> perform same calculation slower, but
#                               more memory friendly
#     """
#     max_label = label_cube.max()
#     if memory_friendly == False:
#         det_sn_mean = measurements.mean(det_sn_cube,
#                                     labels=label_cube,
#                                     index=xrange(1,max_label+1))
#         det_sn_std_dev = measurements.standard_deviation(det_sn_cube,
#                                                      labels=label_cube,
#                                                      index=xrange(1,
#                                                                   max_label+1))
#     else:
#         det_sn_mean = np.empty(max_label)
#         det_sn_std_dev = np.empty(max_label)
#         # instead of calculating for all labels at once,
#         # iterate over labels...
#         for label in xrange(1,max_label+1):
#             det_sn_mean[label - 1] = measurements.mean(det_sn_cube,
#                                                        labels=label_cube,
#                                                        index=label)
#             det_sn_std_dev[label - 1] =\
#                 measurements.standard_deviation(det_sn_cube,
#                                                 labels=label_cube,
#                                                 index=label)
#     return det_sn_mean,det_sn_std_dev

##

# def calc_sn_extrema(det_sn_cube,label_cube,memory_friendly=False):
#     """
#     det_sn_min,det_sn_max,x_snmax,y_snmax,z_snmax,x_snmin,y_snmin,
#     z_snmin = calc_sn_extrema(det_sn_cube,label_cube,memory_friendly=False)

#     Calculation of S/N extrema using scipy.ndimage.measurements.extrema
#     memory_friendly = True -> perform same calculation slower, but
#                               more memory friendly
#     """
#     # ATTENTION - THIS ROUTINE SUCKS BIGTIME!

#     # TODO: measurements.extrema is slow on big cubes--- different solution?!    
#     # already outsourced det_sn_max & det_sn_min if requested alone..
#     # I THINK I CAN DO SOMETHING SIMILAR ... find_objects -> list comprehension magic

#     max_label = label_cube.max()

#     if memory_friendly == False:

#             sn_extrema = measurements.extrema(det_sn_cube,
#                                               labels=label_cube,
#                                               index=xrange(1,max_label + 1))
#             det_sn_min = sn_extrema[0]
#             det_sn_max = sn_extrema[1]
#             det_sn_min_coords = np.asarray(sn_extrema[2],dtype=int)[:,::-1]
#             det_sn_max_coords = np.asarray(sn_extrema[3],dtype=int)[:,::-1]
#             x_sn_max = det_sn_max_coords[:,0]; x_sn_min = det_sn_min_coords[:,0]
#             y_sn_max = det_sn_max_coords[:,1]; y_sn_min = det_sn_min_coords[:,1]
#             z_sn_max = det_sn_max_coords[:,2]; z_sn_min = det_sn_min_coords[:,2]

#     else: 
#         # above call measurements.extrema with index=xrange...  call
#         # sucks a lot of memory and is slow .. iterating over
#         # the objects seems to be more memory friendly ...
#         # probably a topic for the scipy mailinglist
#         det_sn_min = np.empty(max_label)
#         det_sn_max = np.empty(max_label)
#         x_sn_max = np.empty(max_label); y_sn_max = np.empty(max_label)
#         z_sn_max = np.empty(max_label)
#         x_sn_min = np.empty(max_label); y_sn_min = np.empty(max_label)
#         z_sn_min = np.empty(max_label)

#         for label in xrange(1,max_label + 1):
#             sn_extremum = measurements.extrema(det_sn_cube,
#                                               labels=label_cube,
#                                               index=label)
#             det_sn_min[label - 1] = sn_extremum[0]
#             det_sn_max[label - 1] = sn_extremum[1]
#             det_sn_min_coord = np.asarray(sn_extremum[2],dtype=int)[::-1]
#             det_sn_max_coord = np.asarray(sn_extremum[3],dtype=int)[::-1]            
#             x_sn_max[label - 1]  = det_sn_max_coord[0]
#             y_sn_max[label - 1]  = det_sn_max_coord[1]
#             z_sn_max[label - 1]  = det_sn_max_coord[2]
#             x_sn_min[label - 1]  = det_sn_min_coord[0]
#             y_sn_min[label - 1]  = det_sn_min_coord[1]
#             z_sn_min[label - 1]  = det_sn_min_coord[2]

#     return det_sn_min,det_sn_max,x_sn_max,y_sn_max,z_sn_max,\
#         x_sn_min,y_sn_min,z_sn_min
