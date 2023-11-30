from scipy.ndimage import measurements
from lsd_cat_lib import *
import line_em_funcs as lef
import numpy as np

def nb_zmin_zmax_i(sn_cube,x_peak_sn_i,y_peak_sn_i,z_peak_sn_i,thresh_ana,ws=20.):
    """
    z_min, z_max above analysis threshold for detection at
    x_peak_sn_i, y_peak_sn_i, z_peak_sn_i
    """
    # subcube_border_cuts from lsd_cat_lib
    xwin_min, xwin_max, ywin_min, ywin_max, zwin_min, zwin_max = \
                    subcube_border_cuts(x_peak_sn_i, y_peak_sn_i, z_peak_sn_i,
                                        sn_cube, ws=ws)

    subcube = sn_cube[zwin_min:zwin_max,
                      ywin_min:ywin_max,
                      xwin_min:xwin_max]

    subcube_logical = subcube > thresh_ana

    subcube_logical_labels, subcube_num_det = \
                                measurements.label(subcube_logical)

    # obj_label is the voxel-cluster above the analysis threshold
    # belonging to the original detection
    obj_label = subcube_logical_labels[z_peak_sn_i - zwin_min,
                                       y_peak_sn_i - ywin_min,
                                       x_peak_sn_i - xwin_min]
    
    # get slice for obj_label
    slices = measurements.find_objects(subcube_logical_labels,
                                       max_label=obj_label)
    obj_label_slice = slices[-1]

    z_min_i = int(obj_label_slice[0].start) + zwin_min
    z_max_i = int(obj_label_slice[0].stop) + zwin_min

    return z_min_i, z_max_i


def nb_zmin_zmax(sn_cube,x_peak_sn,y_peak_sn,z_peak_sn,thresh_ana,ws=20.):
    """
    as nb_zmin_zmax_i - but for an ensemble of peak coordinates
    """
    z_min = []
    z_max = []
    for x_peak_sn_i,y_peak_sn_i,z_peak_sn_i in zip(x_peak_sn,
                                                   y_peak_sn,
                                                   z_peak_sn):
        z_min_i, z_max_i = nb_zmin_zmax_i(sn_cube,
                                          x_peak_sn_i, y_peak_sn_i, z_peak_sn_i,
                                          thresh_ana, ws=ws)
        z_min.append(z_min_i)
        z_max.append(z_max_i)

    return np.asarray(z_min), np.asarray(z_max)
                                          

def ana_thresh_region_at_peak_i(sn_cube,x_peak_sn_i,y_peak_sn_i,z_peak_sn_i,
                                thresh_ana):
    """
    logical_mask_at_peak = ana_thresh_region_at_peak_i(...)
    """
    peak_layer = sn_cube[z_peak_sn_i,:,:]
    peak_layer_logical = peak_layer > thresh_ana
    peak_layer_labels, peak_layer_num_det = \
                                        measurements.label(peak_layer_logical)

    obj_label = peak_layer_labels[y_peak_sn_i,x_peak_sn_i]

    logical_mask_at_peak = peak_layer_labels == obj_label

    return logical_mask_at_peak


def ana_thresh_region_at_peak(sn_cube, x_peak_sn, y_peak_sn, z_peak_sn,
                              thresh_ana):
    """
    """
    logical_masks_at_peak = []  # a list holding all the logical arrays
    for x_peak_sn_i, y_peak_sn_i, z_peak_sn_i in zip(x_peak_sn,
                                                     y_peak_sn,
                                                     z_peak_sn):
        logical_mask_at_peak = ana_thresh_region_at_peak_i(sn_cube,
                                                           x_peak_sn_i,
                                                           y_peak_sn_i,
                                                           z_peak_sn_i,thresh_ana)
        logical_masks_at_peak.append(logical_mask_at_peak)
        
    return logical_masks_at_peak


def z_diff_calc(z_mins, z_maxs, z_peak_sns, z_diff_max=20):
    assert np.all(z_maxs > z_mins)
    z_diff = z_maxs - z_mins
    z_diff_select = z_diff >= z_diff_max
    z_mins[z_diff_select] = z_peak_sns[z_diff_select] - int(z_diff_max/2.)
    z_maxs[z_diff_select] = z_peak_sns[z_diff_select] + int(z_diff_max/2.)
    z_diff_flag = z_diff_select.astype('int')
    return z_mins, z_maxs, z_diff_flag


def gen_nb_image(cube, z_min, z_max):
    return np.sum(cube[z_min:z_max,:,:],axis=0)


def gen_nb_images(cube, z_mins, z_maxs):
    nb_images = [gen_nb_image(cube, z_min_i, z_max_i)
                 for z_min_i, z_max_i in zip(z_mins,z_maxs)]
    return nb_images
        

def gen_masked_flux_filt_nb_images(logical_masks_at_peak, nb_images):
    masked_flux_filt_nb_images = [logical_mask_at_peak * nb_image
                                  for logical_mask_at_peak, nb_image
                                  in zip(logical_masks_at_peak, nb_images)]
    return masked_flux_filt_nb_images


def sn_2dmom_img(masked_flux_filt_nb_image):
    """
    x_1mom,y_1mom,x_2mom,y_2mom,xy_2mom = sn_2dmom_img(masked_flux_filt_nb_image)

    FIRST & SECOND IMAGE MOMENTS OF "ISOPHOTAL" SELECTED REGION
    """
    indices_y, indices_x = np.indices(masked_flux_filt_nb_image.shape)

    x_1mom = np.sum(masked_flux_filt_nb_image * indices_x) / \
             np.sum(masked_flux_filt_nb_image)
    y_1mom = np.sum(masked_flux_filt_nb_image * indices_y) /  \
             np.sum(masked_flux_filt_nb_image)
    x_2mom = np.sum(masked_flux_filt_nb_image * indices_x**2) / \
             np.sum(masked_flux_filt_nb_image) - x_1mom**2
    y_2mom = np.sum(masked_flux_filt_nb_image * indices_y**2) / \
             np.sum(masked_flux_filt_nb_image) - y_1mom**2
    xy_2mom = np.sum(masked_flux_filt_nb_image * indices_x * indices_y) / \
              np.sum(masked_flux_filt_nb_image) - x_1mom * y_1mom

    return x_1mom,y_1mom,x_2mom,y_2mom,xy_2mom


def sn_2dmoms(masked_flux_filt_nb_images):
    """
    """
    x_1moms = [] ;  y_1moms = [] ; x_2moms = [] ;  y_2moms = [] ; xy_2moms = []
    for masked_flux_filt_nb_image in masked_flux_filt_nb_images:
        x_1mom,y_1mom,x_2mom,y_2mom,xy_2mom = \
                                        sn_2dmom_img(masked_flux_filt_nb_image)
        x_1moms.append(x_1mom)
        y_1moms.append(y_1mom)
        x_2moms.append(x_2mom)
        y_2moms.append(y_2mom)
        xy_2moms.append(xy_2mom)

    return np.asarray(x_1moms), np.asarray(y_1moms), \
        np.asarray(x_2moms), np.asarray(y_2moms), np.asarray(xy_2moms)


def create_kron_calc_image(image_array,x_1mom,y_1mom,sigma_iso,R_max=6.):
    """
    Create the arrays from which the kron radius will be caluclated.
    These arrays are created by selecting a region R_max * sigma_iso
    centered on x_1mom, x_2mom from each image array. In principle the
    nb_images from summation over z_min -> z_max in the continuum
    subtracted fluxcube could be used. But it is more "stable" to use
    the ones from the filtered flux cube...
    """
    r_pix_from_1mom = pixelrad_from_cen(image_array, x_1mom, y_1mom)
    kron_calc_area = r_pix_from_1mom <= R_max * sigma_iso
    kron_calc_image = kron_calc_area * image_array

    return kron_calc_image


def create_kron_calc_images(image_arrays, x_1moms, y_1moms, sigma_isos,
                            R_max=6.):
    kron_calc_images = []
    for image_array, x_1mom, y_1mom, sigma_iso in zip(image_arrays,
                                                     x_1moms, y_1moms,
                                                     sigma_isos):
        kron_calc_image = create_kron_calc_image(image_array,
                                                 x_1mom, y_1mom, sigma_iso)
        kron_calc_images.append(kron_calc_image)

    return kron_calc_images


def calc_rkron(x_1mom, y_1mom, image_array):
    """
    Calculation of Kron radius around x_1mom, y_1mom of image_array.
    """
    r_pix_from_1mom = pixelrad_from_cen(image_array, x_1mom, y_1mom)
    r_kron = np.sum(image_array * r_pix_from_1mom) / np.sum(image_array)
    return r_kron


def calc_rkrons(x_1moms, y_1moms, image_arrays, R_min=3.):
    """
    """
    r_krons_list = [calc_rkron(x_1mom, y_1mom, image_array)
                    for x_1mom, y_1mom, image_array
                    in zip(x_1moms, y_1moms, image_arrays)]
    r_krons = np.asarray(r_krons_list)

    r_kron_flag = r_krons < R_min
    r_krons[r_kron_flag] = R_min

    return r_krons, r_kron_flag.astype('int')


def flux_rkron_circ(mfs_cube, x_1mom, y_1mom, z_min, z_max, r_kron,
                    k_scale=2.5, err_cube=None, delta_lambda=1.25):
    """
    flux(,err) = flux_rkron_circ(mfs_cube, x_1mom, y_1mom, z_min, z_max, r_kron, k_scale,
                 err_cube=None, delta_lambda=1.25):
    (err only if err_cube != None)
    
    """
    flux_mfs_nb_image = gen_nb_image(mfs_cube, z_min, z_max)
    if err_cube is not  None:
        assert mfs_cube.shape == err_cube.shape
        nb_err_image = gen_nb_image(err_cube, z_min, z_max)
        
    r_pix_from_1mom = pixelrad_from_cen(flux_mfs_nb_image,x_1mom, y_1mom)

    r_kron_select = r_pix_from_1mom <= k_scale * r_kron

    flux_sel_image = r_kron_select * flux_mfs_nb_image
    flux = np.sum(flux_sel_image) * delta_lambda
    
    if err_cube is not None:
        err_sel_image = r_kron_select * nb_err_image
        err = np.sqrt(np.sum(err_sel_image)) * delta_lambda

        return flux, err
    else:
        return flux


def fluxes_rkron_circ(mfs_cube, x_1moms, y_1moms, z_mins, z_maxs, r_krons,
                      k_scale=2.5, err_cube=None, delta_lambda=1.25):
    """
    """
    if err_cube is None:
        fluxes = [flux_rkron_circ(mfs_cube, x_1mom, y_1mom, z_min, z_max, r_kron,
                                  k_scale=k_scale, delta_lambda=delta_lambda) 
                  for x_1mom, y_1mom, z_min, z_max, r_kron
                  in zip(x_1moms, y_1moms,z_mins, z_maxs, r_krons)]
        return np.asarray(fluxes)
    else:
        fluxes = []
        errs = []
        for x_1mom, y_1mom, z_min, z_max, r_kron in zip(x_1moms, y_1moms,
                                                        z_mins, z_maxs, r_krons):
            flux, err = flux_rkron_circ(mfs_cube, x_1mom, y_1mom, z_min, z_max, r_kron,
                                        err_cube=err_cube,
                                        k_scale=k_scale, delta_lambda=delta_lambda)
            fluxes.append(flux)
            errs.append(err)

        return np.asarray(fluxes), np.asarray(errs)
        
        
