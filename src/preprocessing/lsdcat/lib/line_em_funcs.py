'''
    File name: line_em_funcs.py
    Author: Christian Herenz

    Convience functions for use in LSDCat
'''

import re
import sys
import time

import numpy as np
from astropy.io import fits


def int_or_str(hdu_str):
    """Helper for argparse.ArgumentParser.add_argument() so that either
    HDU numbers or HDU names can be supplied as commandline arguments.
    """
    if hdu_str.isdigit():
        return int(hdu_str)
    else:
        return hdu_str


def wavel(header, naxis=3, cdel_key='CD3_3'):
    """
    xax = wavel(header,naxis=3,cdel_key='CD3_3')
    header ... a pyfits header object containing naxis,
     crval, crpix & crdel values
               needed for creation of wavelength grid array
    naxis ... the wavelength axis (used to determine CRVALn,
     NAXISn & CRPIXn values)
    cdel_key ... the key word that contains the wavelength
     increment array
    """
    nax = naxis

    naxis = header['NAXIS' + str(nax)]
    crval = header['CRVAL' + str(nax)]
    crpix = header['CRPIX' + str(nax)]
    crdel = header[cdel_key]

    xax = crval + (np.arange(naxis) - (crpix - 1)) * crdel

    return xax


def get_timestring(starttime=None):
    """
    time_string = get_timestring(starttime=None)

    returns empty string, if starttime=None otherwise
    the a string 'xxx.xxx s' where xxx.xxx are the seconds
    that have been ellapsed since starttime (i.e.
    an earlier  time.time() call in your script)
    """
    assert isinstance(starttime, type(time.time()))
    if starttime is None:
        return ''
    else:
        return '(' + str(round(time.time() - starttime, 3)) + 's)'


def autoscale_vmin_vmax(image, scale_factor):
    """
    vmin,vmax = autoscale_vmin_vmax(image,scalefactor)

    Determine vmin,vmax values for matplotlib.imshow (matplotlib.set_ylim)
    such that only scalefactor*100 percent of the values in the image
    (spectrum)

    Image should be either an 1D spectrum or an 2D image
    (other dimensions are also allowed, but do not make
    much sense).
    """
    assert 1 >= scale_factor > 0
    assert isinstance(image, np.ndarray)
    image = np.sort(image.flatten())
    arg_vmax = image.size * scale_factor
    arg_vmin = image.size * (1 - scale_factor)
    vmax = image[int(arg_vmax)]
    vmin = image[int(arg_vmin)]
    return vmin, vmax


def read_hdu(infile, hdunum, nans_to_value=False, nanvalue=0, memmap=True):
    """
    image_data, header = read_hdu(infile,hdunum,nans_to_value=False,nanvalue=0)
    Read in a HDU <hdunum> from FITS-File <infile>, returns the
    data <image_data> and the header <header>.
    Additionally with nans_to_value = True (default: False) one can set
    all NaNs to a specific <nanvalue> (default: 0)
    """
    hdu = fits.open(infile, memmap=memmap)
    header = hdu[hdunum].header
    image_data = hdu[hdunum].data
    hdu.close()
    if nans_to_value is True:
        # change NaNs to nanvalue
        select = np.isnan(image_data)
        image_data[select] = nanvalue
    return image_data, header


def write_primary(data, fheader, outfile):
    """
    Write out primary HDU with <data> and <fheader> into a FITS
     file named <outfile>:
    write_primary(data,fheader,outfile)
    """
    out = fits.PrimaryHDU(data, header=fheader)
    out.writeto(outfile,
                overwrite=True,
                output_verify='silentfix')


def write_fitscube(data,
                   data_header,
                   primary_header,
                   outfile,
                   overwrite=False):
    primary_hdu = fits.PrimaryHDU(data=None, header=primary_header)
    data_hdu = fits.ImageHDU(data=data, header=data_header)
    hdu_list = fits.HDUList([primary_hdu, data_hdu])
    try:
        hdu_list.writeto(outfile, overwrite=overwrite)
    except OSError:
        print("""
ERROR WRITING OUTPUT FITS FILE %s.
Check if file exists (you may use --overwrite for overwriting the
file). If file not exists check if you have write permissions in the
output directory, or if the disk is full. Terminating now!""" % outfile)
        sys.exit(2)


def sexcat_to_reg(sexcat, regfile):
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
    ofile = open(regfile, 'w')

    # Header of regfile:
    ofile.write('# Region file format: DS9 version 4.0 \n')
    ofile.write('\n global dashlist = 5 5 font="helvetica 10 normal"  \n \n')

    for line in ifile:
        if '#' in line:
            continue  # ignore header of ASCII File, ignore comments

        catalog = line.split()
        # ID
        catno = str(catalog[0])
        # Position on the Image
        x = float(catalog[3])
        y = float(catalog[4])
        # Ellipsial Parameters
        kron = float(catalog[2])
        # minor and major axis given in kron radii
        a = kron * float(catalog[5])
        b = kron * float(catalog[6])
        theta = float(catalog[7])
        # Flag, and S/G Classifier
        flag = int(catalog[8])
        star = float(catalog[9])

        # Write the lines of the regfile:
        if flag < 3:
            ofile.write('ellipse(%(x).3f,'
                        ' %(y).3f, %(a).3f,'
                        ' %(b).3f, %(theta).1f) ' % vars())
            if star < 0.5:
                ofile.write('# text = {%(catno).5s} \n' % vars())
            else:
                ofile.write(
                    '# text = {%(catno).5s} color = yellow \n' % vars())
        else:
            ofile.write('ellipse(%(x).3f,'
                        ' %(y).3f, %(a).3f,'
                        ' %(b).3f, %(theta).1f) # dash = 1 ' % vars())
            if star < 0.5:
                ofile.write('text = {%(catno).5s} \n' % vars())
            else:
                ofile.write('text = {%(catno).5s} color = yellow \n' % vars())


def expand_cat(SX_input):
    """
    For the GALFIT fitting, it is neccesary to add a new column to the
    SExtractor output catalog. This column is used for bookkeeping purposes.
    """
    fin = open(SX_input + '.orig', 'r')
    fout = open(SX_input, 'w')
    fout.write('#   0 FITFLAG         If 1 the program will fit this object\n')
    for line in fin:
        if line[0] == '#':
            fout.write(line)
        else:
            fout.write('1' + line)

    fin.close()
    fout.close()


def extract_spectra_circ(cube, x, y, r, extr_func=np.mean):
    """
    spectrum = extract_spectra(fits, x, y, extr_func=np.mean):

    Extracts spectra on datacubes (given in <fits>,<hdu>)
     in circular apertures (well,
    only inside of points which are inside the circle...)
     of radius <r> at postion <x>,<y>.
    Output is 1D-array <spectrum>, with each element
    containing the mean of each spectral-channel in the aperture.

    extr_func controls which method is being used to calculate the
     output spectrum.
    (e.g. could be np.mean [default], np.std, or np.sum,
     or something like this).
    """

    # cube,data = read_hdu(infile,hdu)
    xmax = cube.shape[2]
    ymax = cube.shape[1]
    zmax = cube.shape[0]

    # test for bogus coordinates or no cube
    if x - 1 > xmax or y - 1 > ymax or len(cube.shape) != 3:
        sys.stdout.write('Error!')
        sys.exit(2)
    # TODO: test for radius, that is out of bounds
    # TODO: aperture photometry
    # make boolean array with spaxels to extract = true:
    xcen = x - 1
    ycen = y - 1
    ygrid, xgrid = np.indices((ymax, xmax), float)
    dy = ygrid - ycen
    dx = xgrid - xcen
    radius = np.sqrt(dy ** 2 + dx ** 2)
    select = radius <= r

    # create an array with the fluxes at each spectral "bin"
    spectrum = np.zeros(zmax, dtype=cube.dtype)
    for i in range(zmax):
        spectrum[i] = extr_func(cube[i, :, :][select])

    return spectrum


def hierarch_multi_line(header, keyword, value, comment):
    field_len = 80 - len(keyword) - 4  # assumes keyword inc. "HIERARCH "
    if len(value) <= field_len:
        header[keyword] = (value, comment)
        return header
    else:
        header[keyword] = value[:field_len]
        remaining_val = value[field_len:]
        i = 0
        while len(remaining_val) >= field_len - 1:
            i += 1
            header[keyword + str(i)] = remaining_val[:field_len - 1]
            remaining_val = remaining_val[field_len - 1:]
        header[keyword + str(i + 1)] = (remaining_val, comment)
        return header


def ccs_header(ccs_routine, version, input_filename,
               data_hdu, stat_hdu, pix_scale, mask_filename, mask_hdu,
               filter_name, trunc_constant, pc, lambda_0, bc,
               cube_header):
    cube_header['HIERARCH LSD CCS'] = (
        ccs_routine,
        'spatial cross-correlation (CCS) routine')
    cube_header['HIERARCH LSD CCSV'] = (version,
                                        'CCS version')

    cube_header = hierarch_multi_line(
        cube_header, 'HIERARCH LSD CCSIN', input_filename,
        'CCS input filename')
    cube_header['HIERARCH LSD CCSINS'] = (
        data_hdu,
        'CCS input data HDU - 0-indexed')
    cube_header['HIERARCH LSD CCSINN'] = (
        stat_hdu,
        'CCS input variance HDU - 0-indexed')
    cube_header['HIERARCH LSD CCSPXSC'] = (
        pix_scale,
        'assumed pix scale arcsec/pix in datacube')
    cube_header = hierarch_multi_line(
        cube_header, 'HIERARCH LSD CCMSK', mask_filename,
        'SSC mask filename')
    cube_header['HIERARCH LSD CCSMSKH'] = (
        mask_hdu,
        'CCS mask filename HDU')
    cube_header['HIERARCH LSD CCSFILT'] = (
        filter_name,
        'CCS filter funtion - i.e. Moffat or Gaussian')
    cube_header['HIERARCH LSD CCSFILTT'] = (
        trunc_constant,
        'CCS filter truncated at CCSFILTT x FWHM')
    cube_header['HIERARCH LSD CCSPLY'] = (
        True, 'p(l) = sum_n p_n (l-l0)^n')
    cube_header['HIERARCH LSD CCSPLYL0'] = (lambda_0, 'l0')

    for pi in enumerate(pc[::-1]):
        if pi[0] == 0:
            unit = '[arcsec]'
        elif pi[0] == 1:
            unit = '[arcsec/AA]'
        else:
            unit = '[arcsec/AA**' + str(pi[0]) + ']'

        cube_header['HIERARCH LSD CCSPLYP' + str(pi[0])] = \
            (pi[1], ' '.join(['p_' + str(pi[0]), unit]))

    if filter_name == 'Moffat':
        cube_header['HIERARCH LSD CCSBETA'] = (True,
                                               'beta(l) = sum_n b_n (l-l0)^n')
        for bi in enumerate(bc[::-1]):
            if bi[0] == 0:
                unit = '[1]'
            elif bi[0] == 1:
                unit = '[1/AA]'
            else:
                unit = '[1/AA**' + str(bi[0]) + ']'

            cube_header['HIERARCH LSD CCSB' + str(bi[0])] = \
                (bi[1], ' '.join(['b_' + str(bi[0]), unit]))

    return cube_header


def ccl_header(ccl_routine, version, inputfile, data_hdu,
               noise_hdu, velocity, header, varspec=None):
    header['HIERARCH LSD CCL'] = (ccl_routine,
                                  'spectral cross-correlation (CCL) routine')
    header['HIERARCH LSD CCLV'] = (version, 'CCL version')
    header['HIERARCH LSD CCLINS'] = (data_hdu, 'CCL input data HDU')
    header['HIERARCH LSD CCLINN'] = (noise_hdu, 'CCL input variance HDU')
    header['HIERARCH LSD CCLVFWHM'] = (velocity, 'CCL filter FWHM [km/s]')
    if varspec is not None:
        header = hierarch_multi_line(header,
                                     'HIERARCH LSD CCLEAVS', varspec,
                                     'CCL 1d varspec filename')

    header = hierarch_multi_line(header, 'HIERARCH LSD CCLIN', inputfile,
                                 'CCL input filename')
    return header


def mf_header(mfr_routine, version, width, gwidth, fitscube,
              signalHDU, varHDU, header):
    header['HIERARCH LSD MFR'] = (
        mfr_routine, 'median filter subtract routine')
    header['HIERARCH LSD MFRV'] = (version, 'MFR version')
    header['HIERARCH LSD MFRW'] = (width, 'MFR median filter width')
    header['HIERARCH LSD MFRGW'] = (gwidth, 'MFR gauss sigma width')
    header = hierarch_multi_line(
        header, 'HIERARCH LSD MFRIN', fitscube,
        'ENR input FITS file')
    header['HIERARCH LSD MFRINS'] = (
        signalHDU, 'ENRIN flux HDU name/number')
    if varHDU != -1:
        header['HIERARCH LSD MFRINN'] = (
            varHDU, 'ENRIN variance HDU number')

    return header


def search_header(search_routine, version, inputfile, shdu, nhdu,
                  expmap, expmaphdu, thresh, tabvalues, radius,
                  spaxscale, border_dist, zeroidx, segcube, header):
    header['HIERARCH LSD CAT'] = (search_routine, 'LSDCat search routine')
    header['HIERARCH LSD CATV'] = (version, 'LSDCat search version')
    header = hierarch_multi_line(header, 'HIERARCH LSD CATIN', inputfile,
                                 'LSDCat search infile')
    if nhdu is not None:
        header['HIERARCH LSD SHDU'] = (shdu, 'LSDCat search classic SHDU')
        header['HIERARCH LSD NHDU'] = (nhdu, 'LSDCat search classic NHDU')
    else:
        header['HIERARCH LSD SHDU'] = (shdu, 'LSDCat search MF HDU')

    if expmap is not None:
        header = hierarch_multi_line(header, 'HIERARCH LSD EMAP', expmap,
                                     'LSDCat search exposure map')
        header['HIERARCH LSD EMAPHDU'] = (
            expmaphdu, 'LSDCat search exposure map HDU')

    header['HIERARCH LSD THRESH'] = (
        thresh, 'LSDCat search detection threshold')
    header = hierarch_multi_line(
        header, 'HIERARCH LSD TABVALS', tabvalues,
        'LSDCat search tabvalues')
    header['HIERARCH LSD RADIUS'] = (
        radius, 'LSDCat search grouping radius [asec]')
    header['HIERARCH LSD SCALE'] = (
        spaxscale, 'LSDCat search spaxscale [asec/px]')
    if 'BORDER' in tabvalues:
        header['HIERARCH LSD BDIST'] = (
            border_dist, 'LSDCat search borderdist [px]')

    if zeroidx:
        header['HIERARCH LSD ORIGIN'] = (
            0, 'LSDCat search coords are zero-indexed')
    else:
        header['HIERARCH LSD ORIGIN'] = (
            1, 'LSDCat search coords are one-indexed')

    if segcube is not None:
        header = hierarch_multi_line(
            header, 'HIERARCH LSD SEGCUBE', segcube,
            'LSDCat search segmentation cube')

    return header


def copy_hist_header(header, other_header):
    """
    Copy HISTORY cards from other_header into header,
    but remove duplicates.
    """
    try:
        for hist_item in other_header['HISTORY']:
            header['HISTORY'] = hist_item

        # uniquify history items
        hist_list_o = list(header['HISTORY'])
        hist_list = [i for n, i in enumerate(hist_list_o)
                     if i not in hist_list_o[:n]]

        # remove all items and populate with unique items
        header.remove('HISTORY', remove_all=True)
        for hist_item in hist_list:
            header['HISTORY'] = hist_item

    except KeyError:
        # silently ignore that other_header does not contain any
        # history
        pass

    return header


def string_from_multiline_head(header, keystart):
    string = ''
    extra_keys = [key for key in header.keys()
                  if bool(re.search('\\d\\Z', key)) and
                  key.startswith(keystart)]
    for key in [keystart] + extra_keys:
        string += header[key]
    return string
