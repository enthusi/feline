[//]: # (DO NOT EDIT THE FILE README.md IN THE ROOT OF THIS REPOSITORY)
[//]: # (./doc/gendoc/README_tpl.md IS THE TEMPLATE THAT GENERATES THIS FILE)
[//]: # (THE TEMPLATE GETS CONVERTED INTO ./README.md VIA  ./doc/gen_readme/gen_readme.sh)

# LSDCat - *Line Source Detection and Cataloguing Tool* #

LSDCat is a conceptually simple but robust and efficient detection
package for emission lines in wide-field integral-field spectroscopic
datacubes.  The detection utilises a 3D matched-filtering approach for
compact single emission line objects.  Furthermore, the software
measures fluxes and extents of detected lines.  LSDCat is implemented
in Python, with a focus on fast processing of the large data-volumes
of typical wide-field IFS datacubes.

![LSDCat Logo](https://bitbucket.org/Knusper2000/lsdcat/raw/270da7489147414af85cc191e21a2226981fe6f8/doc/lsd_cat.jpg)

**LSDCat2.0 introduced a major change to the algorithm - read the section
"LSDCat2.0 - What's New" below to learn more.***

The manual is work in progress.  While all the basics of LSDcat are
covered, some of the more advanced features are not yet documented.
Feedback, questions, or even pull requests are more than welcome!


**Table of Contents**

[TOC]

## Requirements ##

LSDCat runs on Python 3.  In order to use LSDCat you need the
following 3rd party Python libraries installed on your system.

  * astropy (>= 4) - http://www.astropy.org/
  * NumPy (>= 1.18) - http://www.numpy.org/
  * SciPy (>= 1.4) - http://www.scipy.org/

A recommended python environment for astronomical research is [astroconda](https://astroconda.readthedocs.io/en/latest/).

## Recommended Additional Software ##

LSDCat can detect emission line objects, but (as of yet) not classify
them.  This functionality is achieved by the companion software
**QtClassify** (http://ascl.net/1703.011). 

QtClassify can automatically identify galaxies that have multiple
emission lines detected.  It also offers an easy-to-use graphical
interface that allows for the interactive classification of single
line detections.  For more information, installation-, and usage
instructions visit the QtClassify web page:
https://bitbucket.org/Leviosa/qtclassify


## Installing LSDCat ##

In your `/home/ehere` (or wherever you want to have LSDCat installed) do:

`git clone https://bitbucket.org/Knusper2000/lsdcat.git`

Then set up your `PYTHONPATH` to contain the LSDCat library directory
`lsdcat/lib`.  For example if you have installed LSDCat in `/home/ehere` and
use bash put in your `.bashrc`

`export PYTHONPATH=${PATH}:${HOME}/lsdcat/lib/`

To have the executables available anywhere on your platform, configure
the systems `$PATH` variable accordingly.  Again, assuming you have
installed LSDCat in your `$HOME` and you are using `bash` add to your
`.bashrc`:

`export PATH=${PATH}:${HOME}/lsdcat/`

If you want to use the additional tools that are shipped with LSDCat
(see below), put

`export PATH=${PATH}:${HOME}/lsdcat/tools/`

in your `.bashrc`.

Users of other shells (e.g., `csh` or `tcsh`) have to follow [a
different procedure](https://kb.iu.edu/d/acar). 

## License ##

LSDCat is licensed under a
[three-clause BSD license](http://choosealicense.com/licenses/bsd-3-clause/). For
details see the file `LICENSE` in the LSDCat repository.

## Acknowledging / Citing LSDCat ##

If your research benefits from the use of LSDCat we ask
you to cite our LSDCat paper:

E. C. Herenz & L. Wisotzki 2016, A&A 602, A111

ADS: https://ui.adsabs.harvard.edu/#abs/2017A%26A...602A.111H
DOI: https://doi.org/10.1051/0004-6361/201629507   (open access)

We also have a record in the Astrophysics Source Code Library:
http://ascl.net/1612.002

More recent version of `LSDCat` (≥1.5) differ slightly from the
algorithm presented in the paper.  Polynomials can now be used to
describe the wavelength dependence of β and the polynomials both for β
and FWHM can be of arbitrary degree.  Moreover, the
`/tools/median-filter-cube.py` routine for an approximate continuum
subtraction can deal with the wavelength range (~5800 Å – 5970 Å)
around NaD (589nm) where the notch-filter blocks light to avoid
contamination and saturation by the laser guide stars.  

With LSDCat2.0 a new search algorithm is the default, this will be
described below.  Both `lsd_cc_spatial.py` and `lsd_cc_spectral.py`
can be run with a `--classic` switch, that restores the old behaviour
and ensures backwards compatibility.  Below we describe the changes a
bit more detailed.  If you make use of the new algorithm you are
kindly asked to cite the following publication:

E.C. Herenz 2022, Astronomical Notes Vol. 344 (5), e20220091

ADS: https://ui.adsabs.harvard.edu/abs/2023AN....34420091H/abstract
DOI: https://doi.org/10.1002/asna.20220091 (open access)

## LSDCat 2.0 - What's new? ##

LSDCat 2.0 implements an updated algorithm for the matched filtering.
An overview of the new algorithm and its motivation has been presented
in a poster at the Joint Observatories Kavli Science Forum in Chile
(ESO Santiago, Chile, April 25-29, 2022): [Link to the poster /
doi:10.5281/zenodo.6471629](https://doi.org/10.5281/zenodo.6471629);
hereafter Poster2.0.  A more detailed and thorough description is presented
in the paper:

E.C. Herenz "Revisiting the Emission Line Source Detection Problem in
Integral Field Spectroscopic Data" (Astronomical Notes, Vol. 344 (5),
e20220091](https://doi.org/10.1002/asna.20220091) ).

In essence, the matched filter now directly incorporates the variances
in the filter design; previously variances were propagated according
to the filtering operation.  Affected by this change are the routines
`lsd_cc_spatial.py` and `lsd_cc_spectral.py`. After subsequent
spatial- and spectral filtering with those routines the output
datacube can be directly interpreted as detection significances,
i.e. the division of signal by propagated noise is not applicable
anymore.  The revised ansatz improves the calculation of the
signal-to-noise of emission lines where the sky-background noise
varies strongly and rapidly as a function of wavelength, i.e close to
(or even on) telluric emission lines.  The new algorithm relies on the
assumption that the positional dependence of the noise can be
sufficiently well described just by rescaling with the exposure map.
For MUSE data this assumption is usually valid for observations that
use typical dithering and rotation patterns.  The user now has to
supply a 1D variance spectrum and an exposure map (optional, if
rescaling with the exposure map is desired).

### Change to spatial convolution routine ###

`lsd_cc_spatial.py` does not need require any new parameters for the
improved procedure, since the spatial filtering step in the new
formalism requires only a change in normalisation (right-hand side of
Eq. 4 in Poster2.0).  Only a single output cube is produced.  This
cube then needs to be supplied as input for the spectral filtering
step.


### Change to spectral convolution routine ###

`lsd_cc_spectral.py` now requires a 1D variance spectrum as input
(supplied to the routine via the new `-v` or `--varspec` parameter).
In order to create such a spectrum we provide, as a proof of concept,
a very simple tool (`tools/ext_avg_var_spec.py`) that calculates an
average variance spectrum for a cicular aperture from a variance
datacube.

The positional dependence of the noise is asserted to be solely caused
by the number of exposures that contribute to a spaxel (we assume here
input exposures of the same length).  For this purpose an exposure
count map, `Nₑₓₚ(x,y)`, can be provided to the routine via the new
`-e` or `--expmap` parameter.  The weighting of the output cube
`SN(x,y,z)` is then computed as `SN'(x,y,z) =
sqrt(Nₑₓₚ(x,y)/max(Nₑₓₚ)) · SN(x,y,z)`.  We advise to correct the
exposure map for the effects of the spatial convolution.  To achieve
this we provide reference implementation that converts a exposure
count cube (e.g., from the MUSE DRP or from MPDAFs `CubeList.combine`)
into an appropriate exposure count map:
`tools/fov_map_from_expcube.py`. Example usage:
`fov_map_from_expcube.py DATACUBE_candels-cdfs-20_v1.0.fits --filter
gaussian --func median -pc 1.4`; for more details consult
`fov_map_from_expcube.py --help`.


### General changes regarding the output format ###

`LSDCat` now copies the header of the primary HDU from the
input FITS files into the primary HDU of the output FITS files.
Information about the processing with the `LSDCat` routines is
appended to that header.  The processed data are stored in subsequent
HDUs.  This behaviour follows conventions established by the [MUSE
DRO](https://ascl.net/1610.004) and
[MPDAF](https://ascl.net/1611.003), where the header of the primary
HDU is used to store relevant meta-data related to observations and
processing.

The names of the HDUs that store the processed data have also been
changed:

  * `lsd_cc_spatial.py`: 2D convolved data in `DATA_2DCC` (propagated
     variances in `STAT_2DCC` in `--classic` mode).
  * `lsd_cc_spectral.py`: 3D matched filter output stored in
    `DATA_3DCC` (propagated variances in `STAT_3DCC` in `--classic` mode).


## Selection of papers utilising LSDCat ##

  * The MUSE-Wide survey: A first catalogue of 831 emission line
	galaxies *Herenz E.C.; Urrutia, T.; Wisotzki, L. et al.* (2017),
	A&A 606, A12 https://doi.org/10.1051/0004-6361/201731055
  * The MUSE-Wide survey: detection of a clustering signal from Lyman α emitters in the range 3 < z < 6 
    *Diener, C.; Wisotzki, L.; Schmidt, K.B. et al.* (2017), MNRAS 471, 3186
	https://doi.org/10.1093/mnras/stx1677
  * Modeling 237 Lyman-α spectra of the MUSE-Wide survey
	*Gronke, M.* (2017), A&A 608, A139
	https://doi.org/10.1051/0004-6361/201731791
  * Nearly all the sky is covered by Lyman-α emission around high-redshift galaxies
	*Wisotzki, L.; Bacon, R.; Brinchman, J. et al.* (2018), Nature 562, 229
	https://doi.org/10.1038/s41586-018-0564-6
  * The MUSE-Wide Survey: survey description and first data release
	*Urrutia, T; Wisotzki, L.; Kerutt, J. et al.* (2019), A&A 624, A141
	https://doi.org/10.1051/0004-6361/201834656
  * Deciphering the Lyman α blob 1 with deep MUSE observations 
	*Herenz, E.C.; Hayes, M. and Scalarta, C.* (2020), A&A 642, A55
    https://doi.org/10.1051/0004-6361/202037464 
  * Mapping the Morphology and Kinematics of a Lyα-selected Nebula at
    z = 3.15 with MUSE *Sanderson, K.~N.; Prescott, M.~M.~K.;
    Christensen, L. et al.* (2022), ApJ 923, 252
	https://doi.org/10.3847/1538-4357/ac3077

See the [list of citations to the LSDCat paper at ADS](https://ui.adsabs.harvard.edu/abs/2017A%26A...602A.111H/citations) for more recent citations.

## Contact ##

For bug reports or feature request you can issue a 
[ticket in the bugtracker](https://bitbucket.org/Knusper2000/lsdcat/issues?status=new&status=open).

You can also simply send an email to `edmund.herenz@iucaa.in`.
Questions, comments, and feedback are always welcome!

## FAQ ##

  - **Q:** Do you have a cat?  **A:** Yes, in fact three cats. Their names are Spaxel, Voxel, and Pixel.

## Documentation ##

### Overview ###

The following flowchart illustrates the processing steps of LSDCat
from an input datacube to a catalogue of positions, shape parameters
and fluxes of emission line sources.

![LSDCat Flowchart](./doc/lsd_cat_flow.png)

Each of the processing steps has an associated LSDCat routine:

- Spatial filtering: `lsd_cc_spatial.py`
- Spectral filtering: `lsd_cc_spectral.py`
- Thresholding: `lsd_cat.py`
- Measurements: `lsd_cat_measure.py`

A complete description of the algorithms can be found in the [LSDCat paper](https://doi.org/10.1051/0004-6361/201629507).  Here we will
describe how to use these tools on wide-field IFS data.

Moreover, we also provide some tools for working with IFS datacubes in
the LSDCat context.  These tools are in the folder `./tools/`.  A
brief overview of the functionality of these tools is given at the end
of the documentation in the section "Additional Tools".

### Input data format ###

LSDCat works with IFS datacubes stored as FITS files.  A FITS file
storing a datacube is assumed to contain two header-data units (HDUs),
one HDU for the flux values and another one for the associated
variances.


### Matched filtering ###

We use a 3D matched filtering approach in LSDCat to obtain a
robust detection statistic for isolated emission line sources in
wide-field IFS datacubes.  Matched-filtering transforms the input
datacube by convolving it with a template that matches the expected
signal of an emission line in the datacube.  LSDCat is primarily
designed for the search of faint compact emission line sources.  For
those sources it is a reasonable assumption that their spatial and
spectral properties are independent.  Therefore, we can perform the 3D
convolution as two successive convolutions, one in each spectral
layer and one along the spectral direction for each spaxel.


#### Spatial filtering ####

For the convolution in each spectral LSDCat offers the choice between
utilising a circular Gaussian profile or a Moffat profile.  Both
functions are commonly used as an approximation of the seeing induced
point spread function (PSF) in ground based optical and near-IR
observations.  The parameter used to characterise the PSF is its full
width at half maximum (FWHM).  The PSF FWHM depends on wavelength.  In
LSDCat the wavelength dependency has to be supplied via the
coefficients of a polynomial

![FWHM(lambda) = sum_n a_n (lambda - lambda_0)^n](./doc/formulae_png/fwhm_lambda.png)

Here the unit of the wavelength is Angstrom, and FWHM is in
arcseconds.  The polynomial coefficients are thus in units of
arcseconds/(Angstrom)^n - where n is the order of the coefficient.
You need to specify these coefficients, as well as the zero-point
`lambda0` in order to run LSDCat.  In the LSDCat paper we describe
several ways on how to determine suitable PSF FWHM coefficients for
your datacubes (see also additional literature presented above).  The
Moffat function includes a second parameter, &beta;, that
parameterises the kurtosis of the PSF.  For observations taken under
natural seeing conditions &beta; can be assumed constant, however for
observations taken with adaptive optics this is not the case.  Thus
if you use a Moffat, &beta; can also be a given as a polynomial

![beta(lambda) = sum_n b_b (lambda - lambda_0)^n](./doc/formulae_png/beta_lambda.png)



##### Note on MPDAF FSF Polynomial Convention ####

If you have obtained the PSF Model polynomials with MPDAF ([MUSE
Python Data Analysis Framework](https://mpdaf.readthedocs.io/)), then
note that MPDAF uses a different convention for the polynomials 
([see
here](https://mpdaf.readthedocs.io/en/latest/muse.html#muse-fsf-models)):

![MPDAF convention](./doc/formulae_png/mpdaf_convention.png)

The conversion from the MPDAF b<sub>i</sub>'s to the LSDCat
a_<sub>i</sub>'s can be achieved via

![conversion formula 1](./doc/formulae_png/conversion_1.png) 

where the c<sub>i</sub>'s follow from the MPDAF b<i>'s

![conversion formula 2](./doc/formulae_png/conversion_2.png) 

using **a** and **b** as shorthands for

![conversion formula a](./doc/formulae_png/conversion_a.png) 

![conversion formula b](./doc/formulae_png/conversion_b.png) 

and

![conversion theta](./doc/formulae_png/conversion_theta.png) 

We provide a method in `./tools/calc_lsd_cat_poly.py` that performs
this conversion.  This script can also read the parameters from the
header of a MPDAF processed cube and print out the relevant part with
the coefficients for the `lsd_cc_spatial.py` call.

##### Usage #####

Spatial filtering in LSDCat is performed by the routine `lsd_cc_spatial.py`.

    :::text
    usage: lsd_cc_spatial.py [-h] -i INPUT [-o OUTPUT] [-S SHDU] [-N NHDU] [--std]
                             [--ignorenoise] [-m MASK] [-M MHDU] [-P PIXSCALE]
                             [-t THREADS] [-bc [BC ...]] [-pc [PC ...]]
                             [--lambda0 LAMBDA0] [--gaussian] [-T TRUNCCONSTANT]
                             [--classic] [--notemp] [--overwrite]
    
    Spatial cross correlation of all the layers in the datacube with a wavelength
    dependent FSF kernel (approximated by a Moffat or a Gaussian).
    
    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Name of the input FITS file containing the flux (and
                            variance) datacube.
      -o OUTPUT, --output OUTPUT
                            Name of the output FITS file. The output FITS file
                            will contain 2 HDUs: In HDU 0 the filtered signal is
                            stored and HDU 1 contains the propagated variances.
                            [Default: `cc_spat_+INPUT`, i.e. `cc_spat_` will be
                            appended to the input file name.]
      -S SHDU, --SHDU SHDU  HDU name or number (0-indexed) in the input FITS file
                            containing the flux data. [Default: 1]
      -N NHDU, --NHDU NHDU  HDU number (0-indexed) or name in the input FITS file
                            containing the variance data. [Default: 2]
      --std                 Some noise cubes are std not variance (e.g. in KMOS).
                            If set the input noise cube will be squared to make it
                            variance.
      --ignorenoise         Switch to not propagate the variance. If set the
                            output FITS file will contain only 1 HDU that stores
                            the filtered signal.
      -m MASK, --mask MASK  Name of a FITS file containing a mask. [Default: none]
      -M MHDU, --MHDU MHDU  HDU number (0-indexed) or name of the mask within
                            MASK-file. The mask is supposed to be a 2D array with
                            the same spatial dimensions containing only ones and
                            zeros. Spaxels corresponding to zero-valued pixels in
                            the mask will be set to zero prior and post the
                            spatial convolution operation. [Default: 1]
      -P PIXSCALE, --pixscale PIXSCALE
                            Size of a spaxel in arcseconds. [Default: 0.2]
      -t THREADS, --threads THREADS
                            Number of CPU cores used in parallel operation.
                            [Default: all available cpu cores]
      -bc [BC ...]          List of polynomial coefficients (whitespace seperated)
                            for the beta parameter of the Moffat profile, i.e. -bc
                            b0 b1 b2 .. bn. For m>n pm == 0. (default: 3.5).
      -pc [PC ...]          List of polynomial coefficients (whitespace seperated)
                            for the PSF FWHM, i.e. -pc p0 p1 ... pn. For m>n pm ==
                            0. (default: 0.8)
      --lambda0 LAMBDA0     Zero-point of the polynomials, in Angstrom. [Default:
                            7050 Angstrom]
      --gaussian            Switch to use a Gaussian profile instead of the
                            default Moffat profile as spatial filter profile. The
                            &beta; parameter will be ignored in that case.
      -T TRUNCCONSTANT, --truncconstant TRUNCCONSTANT
                            Parameter controlling the truncation of the filter
                            window: the filter is truncated at T*WIDTH-PARAM -
                            were WIDTH-PARAM = sigma for Gaussian and FWHM for
                            Moffat. [Default: 8]
      --classic             LSDCat 1.0 mode: Cross-correlate data with filter and
                            propagate the variance; cross-correlated data and
                            propagated varaince will be written to disk.
      --notemp              In --classic mode, do not write out temporary file to
                            free memory prior to propagating the variances. This
                            can speed up things if you have a lot of memory.
      --overwrite           Overwrite output files. Use with caution!
    
    Note: Polynomials for the PSF FWHM and the beta parameter are in the form
    p(lam) = sum_n p_n (lam - lam0)^n, where lambda is in Angstrom and lam0 is
    specified via --lambda0. For the PSF FWHM the p_n's are defined via -pc and
    are assumed to be in units of arcsec/(Angstrom)^n - where n is the order of
    the coefficient. For the dimensionless beta parameter the p_n's are defined
    alike via -bc in units of 1/(Angstrom)^n.

##### Example usage #####

In the following example we want to apply the spatial filtering in a
FITS file `datacube.fits`.  We have determined that the wavelength
dependency of the FWHM can be modelled by a polynomial with
`p0=0.836` arcsec, `p1=-4.4295e-3` arcsec/Angstrom at `lambda0=7050`
Angstrom.  Furthermore, we want to use a 2D Gaussian as a model for
the PSF.

	:::text
	lsd_cc_spatial.py --input=datacube.fits  \
	--gaussian --lambda0=7050 -pc 0.836 -4.4295e-3 --output=spat_c_datacube.fits`

This command will produce a FITS file `spat_c_datacube.fits` that
contains the filtered data in HDU `DATA_2DCC`.  This cube will be used
as input for `lsd_cc_spectral.py`


#### Spectral filtering ####

In LSDCat we adopt as a spectral template a simple 1D Gaussian, where
the width is parameterised by the FWHM in velocity.  The 1D Gaussian
function is an adequate model for the emission lines of unresolved
distant galaxies where often no spatial disentanglement between
ordered motions and unordered motions is possible. 

Implementation details regarding the algorithm are detailed in the document
`./doc/spec_filt_notes/spec_filt_notes.pdf` that is part of this
repository.

##### Usage #####

Spectral filtering in LSDCat is performed by the routine
`lsd_cc_spectral.py`.

    :::text
    usage: lsd_cc_spectral.py [-h] -i INPUT [-v VARSPEC] [-e EXPMAP] [-F FWHM]
                              [-o OUTPUT] [-S SHDU] [-N NHDU] [-E EHDU]
                              [-t THREADS] [--ignorenoise] [--cunit3 CUNIT3]
                              [--nanfile NANFILE] [--nanhdu NANHDU] [--classic]
                              [--notemp] [--nmax NMAX] [--overwrite]
    
    lsd_cc_spectral.py Wavelength smoothing of all spaxels in the datacube with a
    gaussian kernel. Operation is performed on signal and noise HDU.
    
    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Name of the input FITS file containing the spatially
                            cross-correlated datacube (output from
                            lsd_cc_spatial.py).
      -v VARSPEC, --varspec VARSPEC
                            Name of the FITS file containing the variance
                            spectrum. (this option is ignored in --classic mode).
      -e EXPMAP, --expmap EXPMAP
                            Name of the FITS file that contains the exposure map
                            N_exp(x,y) (default: none). If specified, the matched-
                            filter output SN(x,y,z) will be weighted SN'(x,y,z) =
                            sqrt(N_exp(x,y)/max(N_exp)) · SN(x,y,z), where
                            max(N_exp) will be derived from the exposure map.
                            However, you can also fix this value with the --nmax
                            argument (usefull if it is not suitable to extract the
                            variance spectrum in the region of maximum depth).
                            Note that --exmpap is ignored if --classic mode is
                            used.
      -F FWHM, --FWHM FWHM  Specify the FWHM of the Gaussian line template in
                            km/s. [Default: 300 km/s]
      -o OUTPUT, --output OUTPUT
                            Name of the output FITS file. The output FITS file
                            will contain 2 HDUs: In HDU 0 the filtered signal is
                            stored and HDU 1 contains the propagated variances.
                            [Default: `cc_spec_+INPUT`, i.e. `wavelength_smooth_`
                            will be appended to the input file name.
      -S SHDU, --SHDU SHDU  HDU name or number (0-indexed) in the input FITS file
                            containing the flux data. [Default: DATA_2DCC, i.e.
                            the output HDU from lsd_cc_spatial.py]
      -N NHDU, --NHDU NHDU  In normal mode: HDU name or number in the variance
                            spectrum FITS file in which the 1D variances are
                            stored; default 'VARSPEC'. In classic mode: HDU name
                            or number (0-indexed) in the input FITS file
                            containing the variance data; default 'STAT_2DCC'.
      -E EHDU, --EHDU EHDU  HDU name or number (0-indexed) in the exposure map
                            FITS file where the exposure map is actually stored.
                            (This parameter has no effect when --classic mode is
                            used).
      -t THREADS, --threads THREADS
                            Number of CPU cores used in parallel operation.
                            [Default: all available cpu cores]
      --ignorenoise         Classic mode only: Switch to not propagate the
                            variance. If set the output FITS file will contain
                            only 1 HDU that stores the filtered signal. (Note: In
                            normal mode also one HDU is written to disk, but this
                            is correctly normalised to be interpreted as matched
                            filter output.)
      --cunit3 CUNIT3       Specify wavelength unit ('Angstrom' or 'nm' or 'um').
                            [Default: Value from FITS Header.]
      --nanfile NANFILE     Name of an FITS file that contains a 2D image in`
                            --nanmaskhdu` (see below), that is of the same spatial
                            dimensions as the input cube. Spectra corresponding to
                            NaNs in this image will be ignored in the filtering.
                            [Default: None]
      --nanhdu NANHDU       Number of HDU (0-indexed) or name of FITS file
                            specified in --namask, where the 2D image is stored.
                            [Default: 4]
      --classic             LSDCat 1.0 mode: Cross-correlate data with filter and
                            propagate the variance; cross-correlated data and
                            propagated varaince will be written to disk. This
                            requires a variance cube.
      --notemp              In --classic mode, do not write out temporary file to
                            free memory prior to propagating the variances. This
                            can speed up things if you have a lot of memory.
      --nmax NMAX           Weigh the resulting SN(x,y,z) cube via SN'(x,y,z) =
                            sqrt(N_exp(x,y)/nmax)) · SN(x,y,z), where N_exp(x,y)
                            are the values of the exposure map provided with the
                            --expmap parameter. If not set we use SN'(x,y,z) =
                            sqrt(N_exp(x,y)/max(N_exp)) · SN(x,y,z) for the
                            weighting (see --expmap).
      --overwrite           Overwrite output files. Use with caution!

In the context of MUSE IFS data the `--nanfile` option is especially
useful if the datacubes contain a pointing that was observed with a
position angle (PA) significantly different from 0 deg or 90 deg.
This is because the MUSE pipeline samples each observation onto a
rectangular grid where the spatial axis runs from south to north and
where the spectral axis runs from west to east.  Hence, when the PA is
45 deg 50 % of the spaxels within the FITS file will be empty.  The
NaN-mask now allows to ignore these spaxels in the spectral filtering.



##### Example Usage #####

TBD

### Emission line source detection ###

LSDCat detects emission lines by thresholding in the S/N cube which
results from dividing the matched-filtered signal by the propagated
variances. Because of the matched-filtering, the values in the S/N
cube translate into a probability of rejecting the null-hypothesis
that no emission line is present at a given position in the datacube.
This is commonly referred to as the detection significance of a
source.  However, in a strict mathematical sense this direct
translation is only valid for sources that are exactly described by
the matched-filtering template.  Nevertheless, the matched-filtering
performed above always reduces high-frequency noise, while enhancing
sources that are similar to the matched-filter template.  In the
LSDCat paper we quantified the loss of S/N as function of
source-filter mismatch for PSF-like Gaussian emission lines.  There we
explained, that mismatches of the order of 20% between template and
signal result basically in an insignificant reduction of S/N.

The principal input parameter for the emission line source detection
is the detection threshold (`THRESH`).  The above mentioned relation
between threshold and null-hypothesis rejection probability is only
valid if the input variance datacube contains a realistic estimate of
the true noise data.  We recommend, that the detection threshold
should be chosen as the point of diminishing returns after a visual
check of the S/N cube and the distribution of values within it.  A
detection threshold lower than this point will produce a large
increase in spurious detections with only a small compensatory
increase of genuine emission lines.

Thresholding and the construction of the source catalogue is performed
by the routine `lsd_cat_search.py`.  This routine collects all 3D
clusters of neighbouring voxels above the detection threshold.  For
each of these clusters the coordinates of the S/N-peak
(`X_PEAK_SN,Y_PEAK_SN, Z_PEAK_SN`), its value (`DET_SN_MAX`), and the
number of voxels above the detection threshold (`NPIX`) are stored.
In this catalogue each entry also gets assigned a unique identifier (a
so called running ID) `I`.  Moreover, LSDCat can also assign to each
entry an integer object identifier `ID`: multiple detections at a
similar spatial position (within a small search radius, see
description of the `--radius` parameter below) get assigned the same
object identifier.  However, it needs to be checked afterwards whether
these spatial superpositions having the same object identifier are
real objects or two emission line objects at different redshifts.  The
resulting catalogue table is written as a FITS and ASCII table file to
disk.  All pixel coordinates in this output catalogue are 0-indexed.


#### Usage ####

    :::text
    usage: lsd_cat_search.py [-h] -i INPUT [-S SHDU] [-N NHDU] [-e EXPMAP]
                             [-EHDU EXPMAPHDU] [-t THRESH] [-c CATALOGUE]
                             [--tabvalues TABVALUES] [-r RADIUS]
                             [--spaxscale SPAXSCALE] [--borderdist BORDERDIST]
                             [--zeroidx] [--segcube SEGCUBE] [--overwrite]
    
    lsd_cat.py - Create a catalog of sources above a detection
                 threshold in 3D datacube containing signal and
                 variance (presumably a template cross-correlated
                 flux & variance cube from which detection signifcances
                 can be calculated).
    
    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Name of the input FITS file containing either
                            detection significance or matched filtered data and
                            propagated variances.
      -S SHDU, --SHDU SHDU  HDU number (0-indexed) or name of input S/N or
                            filtered cube FITS file that contains the detection
                            significances or matched filtered data. If no `--NHDU`
                            (see below) is supplied, we assume that the input cube
                            in this HDU is S/N, otherwise we assume that it
                            contains the matched filtered data. [Default: 0]
      -N NHDU, --NHDU NHDU  HDU number (0-indexed) or name of filtered cube FITS
                            file that contains the propagated variances of the
                            matched filtering operation. [Default: not set - in
                            this case the datacube in `--SHDU` is interpreted as
                            S/N]
      -e EXPMAP, --expmap EXPMAP
                            FITS file containing a 2D array where the number of
                            exposed voxels are stored for each spaxel. The tool
                            `fov_map_from_expcube.py` can create such a map. This
                            exposure map is required if the output parameter
                            `BORDER` is demanded in the `--tabvalues` option (see
                            below). [Default: None]
      -EHDU EXPMAPHDU, --expmaphdu EXPMAPHDU
                            HDU of exposure map in --expmap. (Default: 'EXP')
      -t THRESH, --thresh THRESH
                            Detection threshold [Default: 8.0]
      -c CATALOGUE, --catalogue CATALOGUE
                            Filename of the output catalogue. Two files will be
                            written on disc, an ASCII catalogue with the specified
                            filename, and an FITS table with `.fits` being
                            appended to this filename. [Default:
                            `catalog_+INPUT+.cat`, where `INPUT` is the name of
                            the input FITS file.]
      --tabvalues TABVALUES
                            Comma-separated list of columns to be written in
                            output catalogue. See below for a list of supported
                            values. [Default:
                            `I,ID,X_PEAK_SN,Y_PEAK_SN,Z_PEAK_SN,DETSN_MAX`]
      -r RADIUS, --radius RADIUS
                            Grouping radius in arcsec. Detections at similar
                            spatial positions within this search radius get
                            assigned the same `ID`. [Default: 0.8 arcsec]
      --spaxscale SPAXSCALE
                            Plate scale (spatial extent of a pixel) in arcsec.
                            [Default: 0.2 arcsec]
      --borderdist BORDERDIST
                            Flag detection in catalogue (column `BORDER`) if it is
                            less than `BORDERDIST` pixels near the field of view
                            border. Only has an effect if the field `BORDER` is
                            requested in `--tabvalues` and if an `--expmap` is
                            supplied. [Default: 10]
      --zeroidx             Switch to zero-indexed output voxel coordinates. (was
                            the default in LSDCat versions < 2).
      --segcube SEGCUBE     Filename of segmentation cube. (Default: None, i.e. no
                            segmenation cube will be written to disk)
      --overwrite           Overwrite already existing output file! USE WITH
                            CAUTION AS THIS MAY OVERWRITE YOUR RESULTS!
    
    Supported values for --tabvalues option:
    I - Running ID 
    ID - Grouped object ID of an entry in the catalog.
    X_PEAK_SN,Y_PEAK_SN, Z_PEAK_SN - Postion of maximum SN.
    LAMBDA_PEAK_SN - same as LAMBDA_0, but for Z_PEAK_SN.
    RA_PEAK_SN, DEC_PEAK_SN - RA & DEC for X_PEAK_SN, Y_PEAK_SN.
    
    NPIX - Number of voxels above detection threshold per source.
    DETSN_MAX - Formal Detection significance of the object
                (i.e. maximum S/N value of the voxels in the obj.).
    BORDER - 1 if objects near the field of view border, 0 otherwise.
            (Note: Distance to border is set via --borderdist parameter 
            [default: 10 px])
            (Note: Distances are calculated wrt. to X_PEAK_SN & Y_PEAK_SN)
    
    NOTE: 
    By default all pixel coordinates in the output catalog are 1-indexed
    (this behaviour can be changed with the --zeroidx switch)
                                    

#### Example usage ####

TBD

### Source parameterisation / measurements ###

LSDCat provides a set of basic parameters for each detection. The
parameters are chosen to be robust and independent from a specific
scientific application.  A detailed description of the available
parameters is given in the LSDCat paper.  For more complex
measurements, involving e.g. fitting of the sources flux
distributions, the LSDCat measurement capability can serve as a
starting point.

Source parameterisation is performed by the routine
`lsd_cat_measure.py`. As input this routine requires the output
catalogue from the detection routine, the matched filtered data
incl. propagated variances, the original data.

The main input parameter influencing the behaviour of the source
parameterisation routine is the analysis threshold (`THRESHANA`).
This additional threshold must be smaller or equal than the detection
threshold.  The role of the analysis threshold in the calculation of
the various parameters is explained in detail in the LSDCat paper.
There we also give guidelines for choosing its value.


#### Usage ####

TBD

#### Example usage ####

TBD

## Advanced Usage ##

In the examples above we presented only the most simple way of using
LSDCat.  Here we now provide an introduction by examples into the more
advanced LSDCat features.

### Subtraction of continuum objects prior matched filtering ###

It is strongly recommended to subtract objects that have detectable
continuum signal within the datacube.  One possibility to remove
continuum signal is to create a datacube median-filtered in spectral
direction and to remove this median-filtered version from the original
datacube.  For this operation we include the tool
`median-filter-cube.py` in the `./tools/` folder.

The following image shows a white-light image of a 1h MUSE cube before
(left) and after (right) the application of `median-filter-cube.py`
with the width parameter being set to `-W 151`.

![Median Filter Example](./doc/median_filter_example.png)


### Utilising a spatial mask in the matched filtering process ###

In some fields it might be benefical to mask out certain regions of
the datacube for better results with LSDCat.  For example, bright
stars or quasars are not well subtracted from the datacube with the
`median-filter-cube.py` tool.  Moreover, in some cases systematic
noise residuals might be present near the borders of the datacube.  To
overcome this issues, a mask can be used in the matched filtering
process.

As an example we show here a white-light image of a MUSE datacube with
2 very bright continuum objects, and some systematic residuals near
the borders. Using the
software [SAO ds9](http://ds9.si.edu/site/Home.html) regions were
drawn that should be excluded in the matched filtering process.


![ds9 region mask example](./doc/ds9.png)


Using the python
library [pyregion](https://github.com/astropy/pyregion) these regions
can be converted into a binary pixel mask.  If the region is saved as
`mask_region.reg` and your datacube (here called
`cube_with_wl_image_in_hdu4.fits`) contains a white-light image, then
the conversion can be done as follows:

	from astropy.io import fits
	import pyregion

	wl_image = fits.getdata('cube_with_wl_image_in_hdu4.fits', 4)
	wl_header = fits.getheader('cube_with_wl_image_in_hdu4.fits', 4)

	mask_regions = pyregion.open('mask_region.reg').as_imagecoord(wl_header)

	mask = mask_regions.get_mask(shape=wl_image.shape)

	mask = ~mask

	fits.writeto('mask.fits', mask.astype('int'))

![Mask image example](doc/mask_image.png)

This mask can now be utilised in the matched filtering with
`lsd_cc_spatial.py`, e.g.

	lsd_cc_spatial.py -i median_filtered_cube_with_wl_image_in_hdu4.fit.fits -m mask.fits


## Additional tools ##

The following set of additional tools are shipped with LSDCat in the
`./tools/` sub-folder.  These scripts provide convenience functions to
pre- or post-process datacube in ways that are usefull in the context
of emission line detection with LSDCat.

- `median-filter-cube.py`: Subtract an in spectral direction
  median-filtered version of the datacube.  Can be used to remove
  sources that have significant detectable continuum signal within the
  datacube.
- `fov_map_from_expcube.py`: Creates from an exposure map datacube an
  exposure map image.  An exposure map datacube contains in every
  voxel the number of exposures that went into this voxel, while an
  exposure map image contains the number of exposure for each spatial
  pixel.  Such a map can be used, e.g., to identify regions without
  any exposures (e.g. field borders). 
- `calc_lsd_cat_poly.py`: Convert FSF coefficients from MPDAF to LSDCat convention.
- `s2n-cube.py`: Create a signal to noise datacube from a FITS file
  containing a signal and a noise HDU.


### Subtraction of running median along spectral axis ###

`median-filter-cube.py` 

    :::text
    usage: median-filter-cube.py [-h] [-S SIGNALHDU] [-V VARHDU] [-o OUTPUT]
                                 [-W WIDTH] [-G GWIDTH] [--passes PASSES] [--ao]
                                 [--notch_start NOTCH_START]
                                 [--notch_end NOTCH_END]
                                 [--notch_bound NOTCH_BOUND] [-t NUM_CPU]
                                 [--savefiltered] [--memmap]
                                 fitscube
    
    Subtract an in spectral direction median-filtered version of the datacube. Can
    be used to remove sources that have significant detectable continuum signal
    within the datacube.
    
    positional arguments:
      fitscube              The FITS File that contains the flux- and variances
                            datacube
    
    options:
      -h, --help            show this help message and exit
      -S SIGNALHDU, --signalHDU SIGNALHDU
                            HDU name or number (0-indexed) containing the flux in
                            fitscube (default: 'DATA').
      -V VARHDU, --varHDU VARHDU
                            HDU name or number (0-indexed) containing the
                            variances in fitscube (default: 'STAT'). Variances
                            will just be copied over to the output cube. Use
                            --varHDU=-1 if this is not desired or if there is no
                            variance cube.
      -o OUTPUT, --output OUTPUT
                            Name of the output FITS file. (default:
                            medfil_W+<width>+<fitscube>)
      -W WIDTH, --width WIDTH
                            The median filter width in spectral pixels. If an even
                            number is given, the actual filter value will be this
                            number +1. (Default: 151)
      -G GWIDTH, --gwidth GWIDTH
                            The sigma of a Gaussian kernel for in spectral pixels.
                            If not 0, the median filtered spaxels will be smoothed
                            with a Gaussian kernel prior to subtraction. This
                            prevents high-frequency leaks. (Default: 0, i.e. no
                            Gaussian post filtering, recommended range: 5 - 10.)
      --passes PASSES       Number of filtering iterations. Subsequent iterations
                            (i.e. median filtering the median filtered data again
                            and again - 2 or 3 times) can sometimes lead to better
                            results (Default: 1).
      --ao                  Fill notch-filter blocked region with linear
                            interpolated values. This prevents an debiases the
                            median filter subtracted datacube close to the edges
                            of the notch filter. See also the parameters
                            --notch_start, --notch_finish, and --notch_bound.
      --notch_start NOTCH_START
                            Start of the notch-filter blocked wavelength range (in
                            Angstrom, default: 5803).
      --notch_end NOTCH_END
                            End of the notch-filter blocked wavelength range (in
                            Angstrom, default: 5970).
      --notch_bound NOTCH_BOUND
                            Size of window before- and after notch-filter blocked
                            region to determine base points for the interpolation.
                            The base point is calculated as the median over this
                            window.
      -t NUM_CPU, --num_cpu NUM_CPU
                            Number of CPU cores to be used (default: all
                            available).
      --savefiltered        Also save the filtered data prior to subtraction.
                            Useful for testing purposes.
      --memmap              Use memory mapping. Reduces memory usage, but might
                            also reduce speed of execution.

### Creation of exposure map from exposure cube ###

`fov_map_from_expcube.py`

    :::text
    usage: fov_map_from_expcube.py [-h] [-E EXPHDU] [--func FUNC]
                                   [--filter FILTER] [-bc [BC ...]] [-pc [PC ...]]
                                   [--lambda0 LAMBDA0] [-t THREADS]
                                   [--zerothresh ZEROTHRESH] [--overwrite]
                                   [--pixscale PIXSCALE]
                                   input
    
    Creates from an exposure map datacube an exposure map image. An exposure map
    datacube contains in every voxel the number of exposures that went into this
    voxel, while an exposure map image contains the number of exposure for each
    spatial pixel. Such a map can be used, e.g., to identify regions without any
    exposures (e.g. field borders). Additionally, the exposure cube can be
    spatially smoothed, which is usefull if the output is to be used as an
    exposure map for lsd_cc_spectral.py.
    
    positional arguments:
      input                 Input FITS cube
    
    options:
      -h, --help            show this help message and exit
      -E EXPHDU, --exphdu EXPHDU
                            FITS HDU name or number (0-indexed) that contains the
                            exposure cube.
      --func FUNC           Function used to compute the exposure map image. Can
                            be either `sum` (default), `mean`, `median`, or 'max'.
      --filter FILTER       The filter used for spatial smoothing, can be either
                            'none' (default), 'gauss', or 'moffat'.
      -bc [BC ...]          List of polynomial coefficients (whitespace seperated)
                            for the beta parameter of the Moffat profile, i.e. -bc
                            b0 b1 b2 .. bn. For m>n pm == 0. Default: 3.5. This
                            parameter has no effect if --filter='none' or
                            --filter='gauss'.
      -pc [PC ...]          List of polynomial coefficients (whitespace seperated)
                            for the PSF FWHM, i.e. -pc p0 p1 ... pn. For m>n pm ==
                            0; Default: 0.8. This parameter has no effect if
                            --filter='none'.
      --lambda0 LAMBDA0     Zero-point of the polynomials, in Angstrom. [Default:
                            7050 Angstrom]. This parameter has no effect if
                            --filter='none'.
      -t THREADS, --threads THREADS
                            Number of CPU cores used in parallel operation.
                            [Default: all available cpu cores]
      --zerothresh ZEROTHRESH
                            Numerical artifacts after FFT convolution can be
                            zeroed-out, i.e. values below this threshold in the
                            final exposure map will be set to 0. This is usefull
                            if the exposure map is being used as a weight map.
                            Default: 0.8.
      --overwrite           Overwrite output files. Use with caution!
      --pixscale PIXSCALE   Linear extent of a pixel/spaxel in arcseconds.
                            Default: 0.2.


### Creation of average variance spectrum from variance cube ###

`ext_avg_var_spec.py`

    :::text
    usage: ext_avg_var_spec.py [-h] -i INPUT [-o OUTPUT] [-N NHDU] [-x XCEN]
                               [-y YCEN] [-r RADIUS]
    
    ext_avg_var_spec.py Extract average variance spectrum from variance HDU;
    needed as input for lsd_cc_spectral.py. Make sure that the aperture covers a
    region where maximum number of exposures, Nₘₐₓ, contributes to the datacube.
    This is needed, as the re-scaling with the 1D variance σ² to a variance map
    σ²(x,y) uses σ²(x,y) = σ² · Nₘₐₓ) / Nₑₓₚ(x,y), where Nₑₓₚ(x,y) is the exposure
    count map. Obviously, this assumes that the integration time for each exposure
    is the same.
    
    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Name of input FITS file containin the variance HDU.
      -o OUTPUT, --output OUTPUT
                            Name of the output FITS file containing the variance
                            spectrum. [Default: eavs_+INPUT, i.e. eavs_ will be
                            appended to the input file name]
      -N NHDU, --NHDU NHDU  HDU number (0-indexed) or name in the input FITS file
                            containing the variance data. [Default: 'STAT']
      -x XCEN, --xcen XCEN  X cooordinate (1-indexed) of the central spaxel of the
                            extraction aperture. [Default: Central spaxel along
                            the x-axis of the datacube]
      -y YCEN, --ycen YCEN  Y cooordinate (1-indexed) of the central spaxel of the
                            extraction aperture. [Default: Central spaxel along
                            the y-axis of the datacube]
      -r RADIUS, --radius RADIUS
                            Radius of the circular extraction aperture (in
                            spaxel). [Default: 10]


### MPDAF - LSDCat polynomial conversion ###

`calc_lsd_cat_poly.py`

    :::text
    
    Usage: calc_lsd_cat_poly.py MPDAF_DATACUBE.fits ref_wl (plot)
           MPDAF_DATACUBE.fits is the FITS datacube that contains a header
           with the FSF* keywords
           ref_wl is the LSDCat reference wavelength that is specified with
           --lambda0 in lsd_cc_spatial.py
    Output: -pc p0 p1 p2 ... pn -bc b0 b1 b2 .. bn  --lambda0 lambda0
            a string that can be directly copied to the call of lsd_cc_spatial.py
    (see also the handwritten notes in ../doc/calc_lsd_cat_poly_notes.jpg
    for the formalism that is used)

### s2n-cube.py ##

    :::text
    usage: s2n-cube.py [-h] -i INPUT [-n NANMASK] [--nanmaskhdu NANMASKHDU]
                       [-o OUTPUT] [-S SHDU] [-N NHDU] [--sigma] [--clobber]
                       [--float64]
    
    Create a signal to noise datacube from a FITS file containing a signal and a
    noise HDU.
    
    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Name of input FITS datacube file (mandatory argument).
                            (default: None)
      -n NANMASK, --nanmask NANMASK
                            FITS file containing a 2D image of same spatial
                            dimensions as cubes. NaNs in this image are then a
                            mask to exclude detections outside the field of view.
                            Best used with with whitelight image. (default: None)
      --nanmaskhdu NANMASKHDU
                            HDU number (0-indexed) containing the NaN mask.
                            (default: 4)
      -o OUTPUT, --output OUTPUT
                            Name of output S/N datacube (default: s2n_+<input>).
                            (default: None)
      -S SHDU, --SHDU SHDU  Name or number of HDU (0-indexed) in input FITS
                            datacube containing the signal. (default: 0)
      -N NHDU, --NHDU NHDU  Name or number of HDU (0-indexed) in input FITS
                            datacube containing the noise. (default: 1)
      --sigma               Switch to interpret noise as sigma. (default: False)
      --clobber             Overwrite already existing output (use at your own
                            risk). (default: False)
      --float64             Write out data in 64 bit. (default: False)
    
    (HDU numbering is 0-indexed)
