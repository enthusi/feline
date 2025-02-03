"""
The `detect_objects` module provides functionalities for converting between
world and pixel coordinates, calculating Gaussian values, and detecting objects
in astronomical images.

Functions:
    - world_to_pix: Converts world coordinates (RA, Dec) to pixel coordinates.
    - pix_to_world: Converts pixel coordinates to world coordinates (RA, Dec).
    - print_lines: Returns catalog lines for a given toggle and redshift.
    - sort_catalog: Sorts the catalog based on a specific value.
    - write_to_file: Writes sorted catalog data to a file.

Dependencies:
    - astropy: For WCS transformations.
    - matplotlib: For visualizations.
    - numpy, scipy: For numerical operations.
    - skimage: For object detection.
    - json, os, sys: For file operations.
"""


import json
import os
import sys
import astropy.io.fits
import astropy.wcs
import imageio
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
import scipy.ndimage
import skimage.feature
import struct

import src.project_path_config as project_path_config


try:
    mpl.use("TkAgg")
except ImportError as e:
    mpl.use("Agg")


# MW23 it is ageneral annoyance to convert pixel positions in an
# image to world coordinates (WCS) given as two angles on the sky (ra, dec)
def pix_to_world(coord: astropy.wcs.WCS, pix: tuple) -> tuple:
    """
    Convert pixel coordinates to world coordinates.

    Given a tuple of pixel coordinates (x, y) and an astropy WCS object,
    this function calculates the corresponding world coordinates in
    right ascension (RA) and declination (Dec).

    Args:
        coord (astropy.wcs.WCS): The World Coordinate System object that
            defines the transformation between pixel and world coordinates.
        pix (tuple): A tuple containing the x and y pixel coordinates.

    Returns:
        tuple: A tuple containing the RA and Dec in degrees.
    """
    pixarray = np.array([[pix[0], pix[1], 0]], np.float64)
    world = coord.wcs_pix2world(pixarray, 0)
    ra = world[0][0]
    dec = world[0][1]
    return ra, dec


def print_lines(toggle: int, z: float) -> list:
    """
    This function returns the lines of the catalog for the given
    toggle and redshift.

    Args:
        toggle: Integer value representing the toggle.
        z: Float value representing the redshift.

    Returns:
        lines: List of strings containing the lines of the catalog.
    """
    lines = []
    for k in range(len(atoms["atoms"])):
        # is k in the template?
        if toggle & 0x1 == 0:
            toggle = toggle // 2
            continue

        # ok, we consider this atom/transition
        toggle = toggle // 2
        atom = atoms["atoms"][k]

        for emission in atom:
            pos = emission * (z + 1)
            name = atoms["atom_id"].get(str(emission))
            lines.append("%s (%.1f)," % (name, pos))
    return lines


def sort_catalog(catalog_lines: list) -> list:
    """
    This function sorts the catalog based on the 5th value in
    descending order.

    Args:
        catalog_lines: List of strings containing the catalog data.

    Returns:
        sorted_catalog: List of strings containing sorted catalog data.
    """
    # Separate lines starting with '#' and lines not starting with '#'
    hash_lines = [line for line in
                  catalog_lines if line.startswith('#')]
    non_hash_lines = [line for line in catalog_lines if not
                      line.startswith('#')]

    # Sort the non-hash lines based on the 5th value
    sorted_non_hash_lines = sorted(non_hash_lines,
                                   key=lambda line: -int(line.split()[4]))

    # Combine the sorted and unsorted lines
    sorted_catalog = hash_lines + sorted_non_hash_lines

    return sorted_catalog


def write_to_file(catalog_list: list) -> None:
    """
    This function writes the sorted catalog to a file named
    'sorted_catalog.txt'.

    Args:
        catalog_list: List of strings containing the catalog data.

    Returns:
        None
    """
    with open(f'{project_path_config.DATA_PATH_RUNTIME_FILES}/sorted_catalog.txt', 'w') as file:
        for line1 in catalog_list:
            file.write(line1 + '\n')


def extract_arrays() -> tuple:
    """
    Extracts grid spacing values from the header of a binary file.

    Reads the first 16 bytes of the file, unpacks the Z, X, and Y spacing
    values as floats, and returns them as integers along with the raw header.

    Returns:
        tuple: (dz, xd, yd, header)
            - dz, xd, yd (int): Spacing values from the header.
            - header (bytes): The raw 16-byte header.
    """
    with open(os.path.join(
            project_path_config.DATA_PATH_PROCESSED,
            "raw_reordered_s2ncube.dat"), "rb") as f:
        header = f.read()[:16]

    dz = struct.unpack("f", header[0:4])[0]
    xd = struct.unpack("f", header[4:8])[0]
    yd = struct.unpack("f", header[8:12])[0]

    return  int(dz), int(xd), int(yd), header

def resize_filters(xd: int, yd: int) -> tuple:
    """
    Resizes filter arrays to the specified dimensions.

    Reads a raw binary file, splits it into four arrays, and resizes each to
    the provided (xd, yd) dimensions.

    Returns:
        tuple: (plane, redshift, template, used)
            - plane, redshift, template, used (ndarray): Resized arrays.
    """
    data = np.fromfile(
        os.path.join(project_path_config.DATA_PATH_RUNTIME_FILES,
                     "feline_float32_array.raw"), dtype="float32")
    plane, redshift, template, used = np.split(data, 4)

    plane.resize((xd, yd))
    redshift.resize((xd, yd))
    template.resize((xd, yd))
    used.resize((xd, yd))

    return plane, redshift, template, used


if __name__ == "__main__":
    hdu = astropy.io.fits.open(
        os.path.join(project_path_config.DATA_PATH_PROCESSED,
                     sys.argv[1]))
    coord = astropy.wcs.WCS(hdu[1].header)

    dz, xd, yd, _ = extract_arrays()

    catalog = [f"# Cube dimensions (z,y,x): {dz}, {xd}, {yd}"]

    with open(os.path.join(project_path_config.DATA_PATH_LOOKUP,
                           "atoms.json"), "r") as data:
        atoms = json.load(data)

    mpl.rcParams["savefig.directory"] = "."

    isize = xd * yd
    size = isize

    plane, redshift, template, used = resize_filters(xd, yd)

    # Scale the floating-point values to the range [0, 1]
    plane_scaled = (plane - np.min(plane)) / (np.max(plane) - np.min(plane))

    # Convert the floating-point values to integers in the range [0, 255]
    plane_uint8 = (plane_scaled * 255).astype(np.uint8)

    # Create a new figure with a specific size
    # (you may adjust the size as needed)
    plt.figure(figsize=(10, 10))

    # Display the image
    plt.imshow(plane_uint8, cmap='jet')

    # Remove the axis
    plt.axis('off')

    # Save the image using imageio.imwrite
    imageio.imsave(f"{project_path_config.DATA_PATH_RUNTIME_FILES}/image.png", plane_uint8)

    if os.path.isfile(f"{project_path_config.DATA_PATH_RAW}/imagemask.png"):
        mymask = imageio.v2.imread(f"{project_path_config.DATA_PATH_RAW}/imagemask.png")

        ny, nx = mymask.shape
        for iy in range(ny):
            for ix in range(nx):
                if mymask[iy, ix] == 0xff:
                    plane[iy, ix] = 0

    else:
        pass
    # print "CAUTION! no valid imagemask.png defined yet!"

    data = scipy.ndimage.gaussian_filter(plane, sigma=1)

    # data=plane
    coordinates = skimage.feature.peak_local_max(
        data, min_distance=1, threshold_abs=50,
        exclude_border=10, num_peaks=100)

    width = 0
    bleft = width
    btop = width
    bright = xd - width
    bbottom = yd - width

    hdu_muse = astropy.io.fits.open(
        f"{project_path_config.DATA_PATH_PROCESSED}/image00.fits",
        memmap=False
    )

    data_muse = hdu_muse[1].data
    nan_sel = np.isnan(data_muse)

    xy = []
    ay = []
    ax = []
    for val in coordinates:
        nx, ny = val
        ax.append(nx)
        ay.append(ny)
        xy.append((nx, ny))

    catalog.append("# " + str(len(xy)))

    run_id = 0
    for hit in xy:  # y,x
        y, x = hit
        if x < 1:
            continue
        y = int(y)
        x = int(x)
        u_i = int(used[y, x])
        z_i = redshift[y, x]
        q_i = plane[y, x]
        t_i = int(template[y, x])
        ra, dec = pix_to_world(coord, (x, y))
        catalog.append("%d %d %d %1.6f %d %d %d" %
                       (run_id, int(y), int(x), z_i, q_i, u_i, t_i) +
                       "\t%.6f %.6f" % (ra, dec) + " " +
                       ' '.join(print_lines(t_i, z_i)))
        run_id += 1

    plt.plot(ay, ax, "*", color="#FFFF00", ms=15,
             path_effects=[path_effects.withStroke(
                 linewidth=3, foreground='black')])
    plt.savefig(project_path_config.DATA_PATH_PDF +
                "/9999_9999_9_f9_fig9999.pdf", format='pdf',
                bbox_inches='tight')

    # Close the figure to free up memory
    plt.close()

    catalog = sort_catalog(catalog)
    write_to_file(catalog)
