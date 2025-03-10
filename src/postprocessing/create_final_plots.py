"""
The `create_final_plots` module generates plots and visualizations for
the cube data. It includes functions for coordinate transformations,
calculating impact parameters, fitting spectral templates, and visualizing data cubes.

Functions:
    - scale_params: Converts redshift to plate scale in kpc per arcsec.
    - get_impact: Calculates the impact parameter between a quasar and a galaxy.
    - get_num_lines: Counts the active emission lines based on a model integer.
    - gauss_function: Computes the value of a Gaussian function.
    - galaxy: Generates a synthetic galaxy spectrum using Gaussian profiles for emission lines.
    - find_bottom_of_plot: Sets axis limits for a plot based on the data range.
    - add_ticks_to_plot: Adjusts plot limits and hides tick labels.

Dependencies:
    - astropy: For cosmological calculations and WCS transformations.
    - matplotlib: For plotting and visualizations.
    - mpdaf: For handling data cubes.
    - numpy: For numerical operations.
    - json, os, sys, logging, struct: For file handling and system operations.
"""


import json
import math
import os
import sys
import astropy.cosmology
import astropy.io.fits
import astropy.wcs
import logging
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import mpdaf
import numpy as np
from scipy.optimize import curve_fit
import ref_index
import struct

import src.project_path_config as project_path_config
import src.postprocessing.detect_objects as detect_objects


def scale_params(redshift: float) -> float:
    """
    Compute the plate scale at a given redshift.

    The plate scale is the distance in kiloparsecs (kpc)
    that corresponds to an angular size of one arcminute
    at the given redshift. This function uses the cosmology model
    defined in the global variable `cosmo`.

    Args:
        redshift (float): The redshift at which to compute the plate scale.

    Returns:
        float: The plate scale at the given redshift, in kpc per arcsec.
    """
    ps = cosmo.kpc_proper_per_arcmin(redshift).value / 60.0
    return ps


# that function gives you the plate scale at the redshift of the galaxy
# to convert arcsec to kpc
# then:
def get_impact(QSO_X: float, QSO_Y: float, px: float, py: float,
               z: float) -> float:
    """
    Calculates the impact parameter, which is the projected distance
    between a quasar (QSO) and a galaxy in kiloparsecs (kpc), at the
    redshift of the galaxy.

    Args:
        QSO_X (float): The x-coordinate of the quasar in arcseconds.
        QSO_Y (float): The y-coordinate of the quasar in arcseconds.
        px (float): The x-coordinate of the galaxy in arcseconds.
        py (float): The y-coordinate of the galaxy in arcseconds.
        z (float): The redshift of the galaxy.

    Returns:
        float: The impact parameter in kpc.
    """

    theta = math.sqrt((QSO_X - px) ** 2 + (QSO_Y - py) ** 2) * 0.2
    scale = scale_params(z)
    # print theta,scale, theta*scale, b
    return theta * scale


# MW23 a single model ist just an integer number,
# here the set bit's are essentially counted
def get_num_lines(toggle: int) -> int:
    """
    Counts the number of active lines (set bits) in the binary
    representation of the input integer. This can be used to determine
    the number of emission lines considered in a given model by
    examining the bits set in the model's integer representation.

    Args:
        toggle (int): An integer where each bit represents the presence
                      (1) or absence (0) of a specific emission line in
                      the model.

    Returns:
        int: The count of active emission lines represented by the set bits in the input integer.
    """
    lines = 0
    for k in range(len(atoms)):
        if toggle & 0x1 == 0:
            toggle = toggle // 2
            continue
        toggle = toggle // 2
        atom = atoms["atoms"][k]
        # atoms_found.append(k)
        for emission in atom:
            lines += 1
    return lines


def gauss_function(x: float, a: float, x0: float, sigma: float) -> float:
    """
    Computes the value of a Gaussian function with given parameters.

    Args:
        x (float): The input value at which to evaluate the function.
        a (float): The amplitude of the Gaussian curve.
        x0 (float): The mean (center) of the Gaussian curve.
        sigma (float): The standard deviation of the Gaussian curve.

    Returns:
        float: The value of the Gaussian function at x.
    """
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


# MW23 used for the red fit to ALL lines at once
# this is now a "proper" galaxy model with
# a Gaussian function for each detected emission
def galaxy(w: float, *p: int) -> np.ndarray:
    """
    Generates a model galaxy spectrum by summing Gaussian profiles for each
    emission line. This function takes a wavelength array and parameters
    for each emission line, including the redshift, line width (sigma),
    and amplitude for each line, to construct a synthetic
    galaxy spectrum. The emission lines are modeled as Gaussian functions.

    Args:
        w (float): An array of wavelengths at which to evaluate the
                   galaxy model.
        p (tuple): A tuple containing parameters for the galaxy model.
                   The first element is the redshift (z), the second is
                   the common line width (sigma) for all lines, followed by
                   the amplitude of each emission line. The number of
                   amplitudes provided should match the number of emission
                   lines being modeled.

    Returns:
        np.ndarray: An array representing the synthetic galaxy spectrum at the wavelengths specified by `w`. This array has the same shape as `w`.
    """
    global forfit_t, atoms
    z = p[0]
    #sigma = p[1]
    sigma=2.0
    # print p
    toggle = forfit_t
    flux = np.zeros(len(w))
    i = 0
    for k in range(len(atoms)):
        if toggle & 0x1 == 0:
            toggle = toggle // 2
            continue
        toggle = toggle // 2
        atom = atoms["atoms"][k]
        # atoms_found.append(k)
        for emission in atom:
            vacline = emission
            pos = vacline * (z + 1)
            newpos = ref_index.vac2air(pos / 10.0) * 10.0
            
            amplitude = p[1 + i]

            flux += gauss_function(w, amplitude, newpos, sigma)

            i += 1
    return flux


def find_bottom_of_plot(ax: plt.Axes, crval:float,
                        crmax: float, data: np.ndarray) -> plt.Axes:
    """
    Sets axis limits for a plot based on the data range.

    Parameters:
        ax (matplotlib.Axes): The plot axes.
        crval (float): Start value for the x-axis.
        crmax (float): Max value for the x-axis.
        data (list or array-like): Data to determine y-axis limits.

    Returns:
        matplotlib.Axes: The modified axes.
    """
    # find bottom:
    lowest = min(data)
    if lowest >= 0:
        bottom = -10
    if lowest < 0:
        bottom = lowest * 1.2

    ax.set_xticks(np.arange(crval, crmax, 200))
    ax.set_xlim(crval, crmax)
    ax.set_ylim(bottom, max(data) * 1.2)

    return ax


def add_ticks_to_plot(plt: plt.Axes, aw: float, px_py: float) -> plt.Axes:
    """
    Adjusts the plot by setting limits, hiding tick labels, and
    inverting the y-axis.

    Parameters:
        plt (matplotlib.pyplot.Axes): The plot axes object to modify.
        aw (float): The center value for the x-axis limits.
        px_py (float): The center value for the y-axis limits.

    Returns:
        matplotlib.pyplot.Axes: The modified axes object.
    """
    plt.xlim(aw - 10, aw + 10)
    plt.ylim(px_py - 10, px_py + 10)
    plt.tick_params(axis="both", which="both", right=True,
                    top=True, left=True, bottom=True,
                    labelbottom=False, labelleft=False,
                    labeltop=False, labelright=False)
    plt.gca().invert_yaxis()

    return plt

def fit_template(t,z,f,w):
    """
    fits a series of galaxy model to a series of lines 
    as indicated by the best match template

    Parameters:
        t (int): The galaxy template 
        z (float): The redshift guess from FELINE
        f (float): The flux array from the median filtered data cube
        w (float): The wavelength array from the median filtered data cube

    Returns:
        new_z (float): The fitted redshift for this galaxy model
    """
    global forfit_t, forfit_w
    forfit_t=t
    forfit_w=w
    params = []
    
    params.append(z)
    #params.append(1.0)
    param_bounds_low=[]
    param_bounds_high=[]
    param_bounds_low.append(z-0.002)
    param_bounds_high.append(z+0.002)
    #sigma bounds in pixels
    #param_bounds_low.append(0.9)
    #param_bounds_high.append(10.0)
    
    #but how many actual lines are that?
    lines=get_num_lines(t)
    for i in range(lines):
        amp=20.0
        sig=1.0
        params.append(amp)
        param_bounds_low.append(0) #amp
        param_bounds_high.append(np.inf)
        
    param_bounds=(param_bounds_low,param_bounds_high)
    popt,pcov = curve_fit(galaxy,w,f,p0=params,bounds=param_bounds,max_nfev=1000)
    new_z=popt[0]
    return new_z

if __name__ == "__main__":
    mpl.use("Agg")
    logging.disable(logging.CRITICAL - 2)
    cosmo = astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

    # These are positions in Angstrom (wavelength)
    # where one might expect absorption in Galaxies! Iron,Magnesium etc.
    # those are just marked as vertical bars later to help the identifcation
    expected_galaxy_absorption_positions = [
        2586.65,
        2600.17,
        2796.35,
        3890.1506,
        3934.777,
        3969.588,
        4102.89,
        4305.61,
        4341.68,
        4862.68,
        5176.7,
        5895.6
    ]

    # place an object, e.g. quasar in pixel coordinates to
    # compute the impact parameter
    qso_x = 150
    qso_y = 150
    foundb = -1
    foundz = -1
    foundm = -1
    foundqop = -1
    foundifrom = "not"

    global px, py, z, run_id, forfit_t, forfit_w
    global used, gtemplate, npx, npy, zerr, ra, dec, quality
    forfit_t = 0
    forfit_w = np.zeros(10)
    global use_new_pos
    use_new_pos = 0

    prev_cat = False

    max_lines_shown = 12
    columns = 12
    rows = 4

    dz, xd, yd, header = detect_objects.extract_arrays()

    crval = struct.unpack("f", header[12:16])[0]
    crmax = crval + dz * 1.25


    if len(sys.argv) < 2:
        print("SYNTAX: %s cube.fits catalog.cat [ds9.reg]" % sys.argv[0])
        sys.exit(0)

    # MW23 lazy way to tell which lines are pairs and
    # which are not something like for:
    # entry in atoms use len(entry)
    atomsize = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1]

    # MW23 just for human readable plotting.
    # This should also be part of the config file
    # SAME positions but now it just gives them a name
    with open(os.path.join(project_path_config.DATA_PATH_LOOKUP,
                           "atoms.json"), "r") as data:
        atoms = json.load(data)

    plane, redshift, template, imused = detect_objects.resize_filters(xd, yd)

    # for data cube
    cube = mpdaf.obj.Cube(os.path.join(project_path_config.DATA_PATH_PROCESSED,
                                       sys.argv[4]), ext=1)
    cubestat = (mpdaf.obj.Cube(os.path.join(
        project_path_config.DATA_PATH_PROCESSED, sys.argv[4]), ext=1))
    # cube.info()

    original_cube = mpdaf.obj.Cube(os.path.join(
        project_path_config.DATA_PATH_RAW, sys.argv[1]), ext=1)

    # MW23 general understanding
    # summing up a CUBE along axis 0 means flattening
    # it into a single image (adds up all layers)
    whiteimage = mpdaf.obj.Cube(os.path.join(
        project_path_config.DATA_PATH_RAW, sys.argv[1])).sum(axis=0).data

    fullwhiteimage = mpdaf.obj.Cube(os.path.join(
        project_path_config.DATA_PATH_RAW, sys.argv[1]), ext=1).sum(axis=0)

    s2ncube = mpdaf.obj.Cube(os.path.join(
        project_path_config.DATA_PATH_PROCESSED, sys.argv[2]), ext=1)
    hdu = astropy.io.fits.open(os.path.join(
        project_path_config.DATA_PATH_PROCESSED, sys.argv[2]))
    coord = astropy.wcs.WCS(hdu[1].header)

    dz, dy, dx = cube.shape

    catalog = open(f"{project_path_config.DATA_PATH_RUNTIME_FILES}/{sys.argv[3]}")
    try:
        colors = mpl.colormaps['winter']
    except Exception as e:
        colors = mpl.cm.get_cmap("winter")
    colors._init()

    i = 0
    vp = 4
    objects = []
    ds = 3

    # interval blue and redward of detection (pix)
    w = 15
    try:
        qso_id = sys.argv[5]

    except Exception as e:
        # print("no QSO given")
        sys.exit(1)

    # MW23 this information is vital but we only
    # have ONE position per data cube
    # it gives the exact position of the central quasar.
    # We should have this function
    # as it contributes alot to certain science cases!
    # should be in config gile, just 3 values in the end:
    # ra,dec of quasar
    # redshif of quasar

    for line in catalog:
        if line[0] == "#":
            continue
        # 496 294 25 1.017953 71 2 16 OIIa (7521.1), OIIb (7526.7),

        use_new_pos = 0

        # reset the values from Lorrie"s catalog in case none is found!
        found = False
        foundb = -1
        foundz = -1
        foundm = -1
        foundqop = -1
        foundifrom = "not"
        # if reading in my own catalog
        if True:
            run_id = int((line.split()[0]))
            py = float((line.split()[1]))
            px = float((line.split()[2]))
            z = float(line.split()[3])

            quality = float((line.split()[4]))
            used = int((line.split()[5]))
            gtemplate = int((line.split()[6]))
            ra, dec = detect_objects.pix_to_world(coord, (px, py))
            # SPECIFIC for THIS cube
            border_distance = min(min(px, py), min(dx - px, dy - py))
            if border_distance < 15:
                continue

        toggle = gtemplate
        positions = []

        atoms_found = []
        lines_found = []

        # MW23 here we sum up over the other 2 axis ->
        # we get a 1d Spectrum (each image becomes ONE pixel essentially)
        # but before we do that, we cut out an area
        # around our object (size is +/- dx and dy...)

        raw_flux = cube[:, int(py) - ds:int(py) + ds,
                        int(px) - ds:int(px) + ds].mean(axis=(1, 2))
        raw_data = raw_flux.data
        raw_wave = np.arange(raw_flux.wave.get_crval(),
                             raw_flux.wave.get_crval() +
                             raw_flux.wave.get_step() *
                             raw_flux.wave.shape,
                             raw_flux.wave.get_step())
        #print("redshift coarse: ",z)
        try:
    
            zfit=fit_template(gtemplate,z,raw_data,raw_wave)
            #print 1/0
        except:
            #print("No galaxy model fitted!")
            zfit=z
        z=zfit
        
        #print("redshift refined: ",z)
        for k in range(len(atoms["atoms"])):
            # is k in the template?
            if toggle & 0x1 == 0:
                toggle = toggle // 2
                continue

            # ok, we consider this atom/transition
            toggle = toggle // 2
            atom = atoms["atoms"][k]
            atoms_found.append(k)
            for emission in atom:
                lines_found.append(emission)
                pos = emission * (z + 1)
                positions.append(pos)

        lines_found.sort()

        wavemin = 4780
        wavemax = 9300

        j = 0
        count = len(positions)
        c = 299792.458
        zguess = z

        # MW23 here I start that one big plot #
        plt.figure(figsize=(16, 9))

        ax1 = plt.subplot2grid((rows, columns), (0, 0), colspan=9)
        plt.title("%s id=%04d, x=%.1f,"
                  " y=%.1f, ra=%.6f dec=%.6f z=%.6f" %
                  (qso_id, run_id, px, py, ra, dec, z))

        ax1.text(0.5, -0.2, r"wavelength ($\rm{\AA}$)", ha="center", va="top",
                 transform=ax1.transAxes, fontsize=10)

        ax2 = plt.subplot2grid((rows, columns), (1, 0), colspan=9)
        plt.title("%d used lines, match strength=%d, b=%.1f" % (
            used, quality, get_impact(qso_x, qso_y, px, py, z)))

        ax2.text(0.5, -0.2, r"wavelength ($\rm{\AA}$)", ha="center", va="top",
                 transform=ax2.transAxes, fontsize=10)

        ax2.tick_params(
            axis="both",  # changes apply to the x-axis
            which="both",  # both major and minor ticks are affected
            bottom="on",  # ticks along the bottom edge are off
            top="off",  # ticks along the top edge are off
            labelbottom="off",
            right="off",
            left="on",
            labelleft="on")  # labels along the bottom edge are off

        # plot regions of absorption first
        # MW23 vac2air is important to compute the position
        # once it was observed THROUGH Earth atmosphere
        for absline in expected_galaxy_absorption_positions:
            abs_wav = ref_index.vac2air(absline * (z + 1) / 10.0) * 10.0
            if crval < abs_wav < crmax:
                ax1.axvline(x=abs_wav, color="aquamarine",
                            linestyle="-", linewidth=4.0)

        for absline in expected_galaxy_absorption_positions:
            abs_wav = ref_index.vac2air(absline * (z + 1) / 10.0) * 10.0
            if crval < abs_wav < crmax:
                ax2.axvline(x=abs_wav,
                            color="aquamarine",
                            linestyle="-",
                            linewidth=4.0)

        # plot actual flux spectrum
        spec = cube[:, int(py) - ds:int(py) + ds,
                    int(px) - ds:int(px) + ds].mean(axis=(1, 2))
        original_spec = original_cube[:, int(py) - ds:int(py) + ds,
                                      int(px) - ds:int(px) + ds].mean(
            axis=(1, 2))
        data1 = spec.data
        original_data1 = original_spec.data
        waven = np.arange(spec.wave.get_crval(),
                          spec.wave.get_crval() +
                          spec.wave.get_step() *
                          spec.wave.shape,
                          spec.wave.get_step())
        ax1.step(waven, original_data1, where="mid", color="darkgrey")

        # plot possible positions for emission

        for a in atoms["atoms"]:
            for b in a:
                p = ref_index.vac2air(b * (z + 1) / 10.0) * 10.0
                if crval < p < crmax:
                    ax1.axvline(x=p, color="k", linestyle="--")

        # plot detected emission  above that
        for g in range(len(positions)):
            wave = lines_found[g]

            thision = atoms["atom_id"].get(str(wave))
            wobs = ref_index.vac2air(wave * (z + 1) / 10.0) * 10.0
            ax1.axvline(x=wobs, color="r", linestyle="--")

        # plot possible positions for emission
        for a in atoms["atoms"]:
            for b in a:
                # p=b*(z+1)
                p = ref_index.vac2air(b * (z + 1) / 10.0) * 10.0
                if crval < p < crmax:
                    ax2.axvline(x=p, color="k", linestyle="--")

        # plot detected emission  above that
        for g in range(len(positions)):
            wave = lines_found[g]
            thision = atoms["atom_id"].get(str(wave))
            wobs = ref_index.vac2air(wave * (z + 1) / 10.0) * 10.0
            # gen narrow bands
            f = 4  # narrow band width
            d = 10
            p = int((wobs - crval) / 1.25)
            all_band = cube[p - f:p + f, :, :]
            if g == 0:
                all_ima = all_band.sum(axis=0)
            else:
                all_ima += all_band.sum(axis=0)

            # wobs=ref_index.vac2air(wobs/10.0)*10.0
            ax2.axvline(x=wobs, color="r", linestyle="--")

        waven_high = np.arange(spec.wave.get_crval(),
                               spec.wave.get_crval() +
                               spec.wave.get_step() *
                               spec.wave.shape,
                               spec.wave.get_step() / 10.0)

        ax1.step(waven, data1, where="mid", color="blue")

        ax1 = find_bottom_of_plot(ax1, crval, crmax, data1)

        s2nspec = s2ncube[:, int(py) - ds:int(py) + ds,
                          int(px) - ds:int(px) + ds].mean(axis=(1, 2))
        data2 = s2nspec.data
        ax2.step(waven, data2, where="mid")

        ax2 = find_bottom_of_plot(ax2, crval, crmax, data2)

        lines_found.sort()

        # plot all found lines AND always Ha,Hb
        first = True
       

        for h in range(min(len(positions), max_lines_shown)):
            ax3 = plt.subplot2grid((rows, columns), (2, h))

            plt.title("%s" % (atoms["atom_id"].get(
                str(lines_found[h]))), fontsize=10)
            wave = lines_found[h]
            # remember Oiii ratio

            #  oxy1=
            # thision = atom_id[wave]
            thision = atoms["atom_id"].get(str(wave))
            wobs = ref_index.vac2air(wave * (z + 1) / 10.0) * 10.0

            # vertical red line in middle of x-axis
            ax3.axvline(x=wobs, color="r", linestyle="--")

            ax3.plot(waven, data1,
                     linestyle="-",
                     drawstyle="steps-mid",
                     color="blue")
            ax3.plot(waven, data2,
                     linestyle="-",
                     drawstyle="steps-mid",
                     color="orange")

            # Set y-axis limits based on the blue line (data1)
            
            

            dl = 15.0
            lim_low = max(crval, wobs - dl)
            lim_high = min(wobs + dl, crmax)
            while (lim_high - lim_low) < (2 * dl):
                lim_high += dl / 3.0
            ax3.set_xlim(lim_low, lim_high)
            ax3.set_xticks([wobs])

            values_zoom=[]
            for i in range(len(waven)):
                wavei=waven[i]
                if wavei>lim_low and wavei<=lim_high:
                    values_zoom.append(data1[i])

            
            y_max = max(values_zoom)*2
            y_min = min(values_zoom)-y_max*0.1
            ax3.set_ylim(y_min, y_max)
            
            
            # only plot axis/label for the first window
            if first:
                ax3.tick_params(
                    axis="both",  # changes apply to the x-axis
                    which="both",  # both major and minor ticks are affected
                    bottom="on",  # ticks along the bottom edge are off
                    top="off",  # ticks along the top edge are off
                    labelbottom="on",
                    right="off",
                    left="on",
                    labelleft="on")  # labels along the bottom edge are off

            else:
                ax3.tick_params(
                    axis="both",  # changes apply to the x-axis
                    which="both",  # both major and minor ticks are affected
                    bottom="on",  # ticks along the bottom edge are off
                    top="off",  # ticks along the top edge are off
                    labelbottom="on",
                    right="off",
                    left="off",
                    labelleft="off")  # labels along the bottom edge are off

            first = False

        # Adjust the space between subplots
        plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05,
                            wspace=0.6, hspace=0.7)

        mark = 0

        bigpic = plt.subplot2grid((rows, 12), (0, 9), rowspan=2,
                                  colspan=4)
        bigpic.imshow(plane, vmax=1000, interpolation="none", cmap="jet")
        bigpic.plot(px, py, "*", color="#FFFF00", ms=15,
                    path_effects=[path_effects.withStroke(
                        linewidth=3, foreground='black')])

        aw = 20
        aw = int(min(aw, px, py))

        spic = plt.subplot2grid((rows, columns), (3, 6))

        # the following block is for the 2nd image, but we need the peak first

        wcs1 = fullwhiteimage.wcs
        narrows = mpdaf.obj.Image(data=all_ima.data,
                                  wcs=wcs1)[int(py) -
                                            aw // 2:int(py) +
                                            aw // 2, int(px) -
                                            aw // 2:int(px) +
                                            aw // 2]

        center_area = narrows[aw // 2 - 3:aw // 2 + 3, aw // 2 - 3:aw // 2 + 3]
        center_mean = np.mean(center_area.data)
        center_std = np.std(center_area.data)
        center_value = center_mean + center_std * 2.0

        plt.imshow(narrows.data, interpolation="none",
                   cmap="jet", vmax=center_value)
        plt.tick_params(axis="both", which="both",
                        right=True, top=True, left=True,
                        bottom=True, labelbottom=False,
                        labelleft=False, labeltop=False,
                        labelright=False)
        plt.title("collapsed")

        # whiteimage plot
        whitezoom = whiteimage[int(py) - aw:int(py) +
                               aw, int(px) - aw:int(px) + aw]
        spic = plt.subplot2grid((rows, columns), (3, 7))

        center_area = whitezoom[aw - 4:aw + 4, aw - 4:aw + 4]
        center_mean = np.mean(center_area.data)
        center_std = np.std(center_area.data)
        # center_peak = whitezoom[aw, aw]
        center_value = center_mean + center_std * 4.0

        plt.imshow(whitezoom, interpolation="none",
                   cmap="jet", vmax=center_value)
        plt.tick_params(axis="both", which="both",
                        right=True, top=True,
                        left=True, bottom=True,
                        labelbottom=False,
                        labelleft=False, labeltop=False,
                        labelright=False)
        plt.title("white")

        plt = add_ticks_to_plot(plt, aw, aw)

        full_plane = mpdaf.obj.Image(data=plane,
                                     wcs=wcs1)[int(py) -
                                               aw:int(py) +
                                               aw, int(px) -
                                               aw:int(px) + aw]

        # quality plot
        spic = plt.subplot2grid((rows, columns), (3, 8))
        center_value = full_plane[aw, aw]

        plt.imshow(full_plane.data, interpolation="none",
                   cmap="jet", vmax=1.0 * center_value)
        plt.title("quality")

        plt = add_ticks_to_plot(plt, aw, aw)

        bestgauss = full_plane.gauss_fit(circular=True,
                                         pos_min=[aw - 4, aw - 4],
                                         pos_max=[aw + 4, aw + 10],
                                         unit_center=None, unit_fwhm=None)
        a, b = bestgauss.center
        fwhm = bestgauss.fwhm

        # redshift plot
        spic = plt.subplot2grid((rows, columns), (3, 9))
        testarea = redshift[int(py) - aw:int(py) + aw,
                            int(px) - aw:int(px) + aw]
        plt.xlim(aw - 10, aw + 10)
        plt.ylim(aw - 10, aw + 10)
        plt.gca().invert_yaxis()
        plt.imshow(testarea, interpolation="none", cmap="jet")
        plt.tick_params(axis="both", which="both", right=True,
                        top=True, left=True, bottom=True,
                        labelbottom=False, labelleft=False,
                        labeltop=False, labelright=False)
        plt.title("redshift")

        # no lines plot
        spic = plt.subplot2grid((rows, columns), (3, 10))
        plt.imshow(imused, interpolation="none", cmap="jet")

        plt = add_ticks_to_plot(plt, px, py)

        plt.title("no. lines")

        npx = px - aw + b
        npy = py - aw + a

        file = "%04d_%04d_%d_f%d_fig%04d.pdf" % (quality,
                                                 run_id,
                                                 mark,
                                                 found,
                                                 i)
        plt.savefig(os.path.join(project_path_config.DATA_PATH_PDF,
                                 file), format="pdf")
        #plt.show()
        plt.close()
        i += 1
