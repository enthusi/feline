import numpy as np
from astropy import wcs as astWCS

__version__ = '1.0.2'

def pix_to_radec(x_coords,y_coords,header): # taken from Christian
    """
    ra_coords, dec_coords = pix_to_radec(x_coords,y_coords,header)
    Calcualtes RA's & DEC's for x- & y-coordinates of datacube.
    header ... header of datacube.
    x_coords ... array containing the 
                 (sub-)x-pixelcoordinates (zero indexed)
    y_coords ... array containing the 
                 (sub-)y-pixelcoordinates (zero indexed)
    """
    wcs = astWCS.WCS(header)
    coords = wcs.wcs_pix2world(x_coords, y_coords,
                               np.zeros(len(y_coords)),0)

    ra_coords = coords[0]
    dec_coords = coords[1]

    return ra_coords, dec_coords

def radec_to_pix(ra_coords,dec_coords,header):
    """
    Reverse of pix_to_radec
    """
    wcs = astWCS.WCS(header)
    world_coords = np.asarray(zip(ra_coords,dec_coords))
    coords = wcs.wcs_world2pix(ra_coords, dec_coords,
                               np.zeros(len(ra_coords)),0)

    x_coords = coords[0]  # I know - this is a little bit strange... 
    y_coords = coords[1]


    return x_coords, y_coords

def spec_in_ap(xcen,ycen,data,ext_radius,xgrid,ygrid):
    # create mask in shape of a circle
    # give summed spectra in the aperture
    # xpos and ypos in pixel: coordinates of an object 
    # ext_radius in pixel
    ##### It follows: Christians Code

    # und so extrahiert man dann ein spektrum
    dy = ygrid - ycen
    dx = xgrid - xcen
    radius2 = np.sqrt(dy**2 + dx**2)

    select = radius2 <= float(ext_radius)
    spec_signal = data[:,select].sum(axis=1)
    return spec_signal
