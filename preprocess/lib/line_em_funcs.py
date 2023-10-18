# Convenience functions for project
# "Detection and Analysis of Line-Emitters in MUSE-Datacubes"

import math
import sys
import time
from astropy.io import fits
import numpy as np

def wavel(header,naxis=3,cdel_key='CD3_3'):
   """
   xax = wavel(header,naxis=3,cdel_key='CD3_3')
   header ... a pyfits header object containing naxis, crval, crpix & crdel values
              needed for creation of wavelength grid array
   naxis ... the wavelength axis (used to determine CRVALn, NAXISn & CRPIXn values)
   cdel_key ... the key word that contains the wavelength increment array
   """
   nax = naxis

   naxis = header['NAXIS'+str(nax)]
   crval = header['CRVAL'+str(nax)]
   crpix = header['CRPIX'+str(nax)]
   crdel = header[cdel_key]

   xax = crval + (np.arange(naxis) - (crpix - 1))*crdel

   return xax

def get_timestring(starttime=None):
    """
    time_string = get_timestring(starttime=None)

    returns empty string, if starttime=None otherwise
    the a string 'xxx.xxx s' where xxx.xxx are the seconds
    that have been ellapsed since starttime (i.e.
    an earlier  time.time() call in your script)
    """
    assert type(starttime) == type(time.time())
    if starttime == None:
        return ''
    else:
        return '('+str(round(time.time()-starttime,3))+'s)'


def autoscale_vmin_vmax(image,scale_factor):
    """
    vmin,vmax = autoscale_vmin_vmax(image,scalefactor) 
    
    Determine vmin,vmax values for matplotlib.imshow (matplotlib.set_ylim)
    such that only scalefactor*100 percent of the values in the image
    (spectrum)

    Image should be either an 1D spectrum or an 2D image
    (other dimensions are also allowed, but do not make
    much sense).
    """
    assert scale_factor <= 1 and scale_factor > 0
    assert type(image) == type(np.ndarray(0))
    image = np.sort(image.flatten())
    arg_vmax = image.size * scale_factor
    arg_vmin = image.size * ( 1 - scale_factor )
    vmax = image[int(arg_vmax)]
    vmin = image[int(arg_vmin)]
    return vmin,vmax


def read_hdu(infile,hdunum,nans_to_value=False,nanvalue=0,memmap=True):
    """
    image_data, header = read_hdu(infile,hdunum,nans_to_value=False,nanvalue=0)
    Read in a HDU <hdunum> from FITS-File <infile>, returns the
    data <image_data> and the header <header>.
    Additionally with nans_to_value = True (default: False) one can set
    all NaNs to a specific <nanvalue> (default: 0)
    """
    hdu = fits.open(infile,memmap=memmap)
    header = hdu[hdunum].header
    image_data = hdu[hdunum].data 
    hdu.close()
    if nans_to_value == True:
        # change NaNs to nanvalue
        select = np.isnan(image_data)
        image_data[select] = nanvalue
    return image_data,header


def write_primary(data,fheader,outfile):
    """
    Write out primary HDU with <data> and <fheader> into a FITS file named <outfile>:
    write_primary(data,fheader,outfile)
    """
    out = fits.PrimaryHDU(data,header=fheader)
    out.writeto(outfile, clobber='True',output_verify='silentfix') # clobber='True' -> overwrite existing file!

def sexcat_to_reg(sexcat,regfile):
    # TODO: add parser for automatic collumn detection
    # TODO: maybe use asciitable or sth... this is not good atm
    """
    Read in an SExtractor ASCII catalog <sexcat> and convert it to a ds9-
    region file <regfile>. The parameter file must be set up in a way,
    such that the first 10 collumns in the ASCII are:
    NUMBER, MAG_AUTO, KRON_RADIUS, X_IMAGE, Y_IMAGE, A_IMAGE, B_IMAGE,
    THETA_IMAGE, FLAGS, CLASS_STAR:
    sexcat_to_reg(sexcat,regfile)
    """
    ifile = open(sexcat, 'r')
    ofile = open(regfile,'w')

    #Header of regfile:
    ofile.write('# Region file format: DS9 version 4.0 \n')
    ofile.write('\n global dashlist = 5 5 font="helvetica 10 normal"  \n \n')

    for line in ifile:
        if '#' in line:
            continue # ignore header of ASCII File, ignore comments

        catalog = line.split()
        # ID
        catno = str(catalog[0])
        # Position on the Image
        x=float(catalog[3]); y=float(catalog[4])
        # Ellipsial Parameters
        kron=float(catalog[2])
        # minor and major axis given in kron radii
        a=kron*float(catalog[5]); b=kron*float(catalog[6])
        theta=float(catalog[7])
        # Flag, and S/G Classifier
        flag=int(catalog[8]); star=float(catalog[9])
    
        # Write the lines of the regfile:
        if flag < 3:
            ofile.write('ellipse(%(x).3f, %(y).3f, %(a).3f, %(b).3f, %(theta).1f) ' % vars())
            if star < 0.5:
                ofile.write('# text = {%(catno).5s} \n' % vars())
            else:
                ofile.write('# text = {%(catno).5s} color = yellow \n' % vars())
        else:
            ofile.write('ellipse(%(x).3f, %(y).3f, %(a).3f, %(b).3f, %(theta).1f) # dash = 1 ' % vars())
            if star < 0.5:
                ofile.write('text = {%(catno).5s} \n' % vars())
            else:
                ofile.write('text = {%(catno).5s} color = yellow \n' % vars())


def expand_cat(SX_input):
    """
    For the GALFIT fitting, it is neccesary to add a new column to the 
    SExtractor output catalog. This column is used for bookkeeping purposes.
    """
    fin = open(SX_input+'.orig', 'r')
    fout = open(SX_input,'w')
    fout.write('#   0 FITFLAG         If 1 the program will fit this object\n')
    for line in fin:
        if line[0] == '#':
            fout.write(line)
        else:
            fout.write('1'+line)

    fin.close()
    fout.close()


def extract_spectra_circ(infile,hdu,x,y,r):
    """
    spectrum = extract_spectra(fits,x,y,a,b,r):
    Extracts spectra on datacubes (given in <fits>,<hdu>) in circular apertures (well,
    only inside of points which are inside the circle...) of radius <r> at postion <x>,<y>.
    Output is 1D-array <spectrum>, with each element
    containing the mean of each spectral-channel in the aperture. 
    """

    cube,data = read_hdu(infile,hdu)
    xmax = cube.shape[2]
    ymax = cube.shape[1]
    zmax = cube.shape[0]

    # test for bogus coordinates or no cube
    if x-1 > xmax or y-1 > ymax or len(cube.shape) != 3:
        sys.stdout.write('Error!')
        sys.exit(2)
    # TODO: test for radius, that is out of bounds
    # TODO: aperture photometry

    # make boolean array with spaxels to extract = true:
    xcen=x-1; ycen=y-1
    ygrid,xgrid = np.indices((ymax,xmax),float)
    dy = ygrid - ycen
    dx = xgrid - xcen
    radius = np.sqrt(dy**2 + dx**2)
    select=radius<=r
    
    # create a table with the fluxes at each spectral "bin"
    spectrum = np.zeros(zmax)
    for i in xrange(zmax):
        spectrum[i] = np.mean(cube[i,:,:][select])

    return spectrum

