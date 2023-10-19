# feline
in parent dir:<br>
wget http://martinwendt.de/cube.fits

## filtering of the raw data cube (median + template filter, see LSDcat)

in preprocess directory:<br>
### apply median filter to 'flat' out the data (emission lines remain):
python2 median-filter-cube.py ../cube.fits --signalHDU=1 --varHDU=2 --num_cpu=6 --width=151 --output=med_filt.fits<br>
### filter the data cube with a 'typical line' template in spatial dimension first:
python2 lsd_cc_spatial_mask.py --input=med_filt.fits --SHDU=0 --NHDU=1 --threads=4 --gaussian --lambda0=7050 -p0=0.7 --output=spatial_cc.fits<br>
### filter the data cube with a 'typical line' template in spectral dimension:
python2 lsd_cc_spectral.py --input=spatial_cc.fits --threads=2 --FWHM=250 --SHDU=0 --NHDU=1 --output=spectral_cc.fits<br>
### filter the data cube with a 'typical line' template in spectral dimension:
python2 s2n-cube.py --input=spectral_cc.fits --output=s2n_v250.fits --clobber<br>
### construct a signal-to-noise cube:
copy s2n_v250.fits and med_filt.fits into parent dir<br>

## here the actual tool starts

in parent dir:<br>

### creating an overview image (for optional masking):
python2 createnan.py cube.fits<br>
### change data order from image-by-image to spectrum-by-spectrum (fast cache access):
python2 preprocess_s2n_cube.py s2n_v250.fits none<br>
### build and run the main Feline code (C/OpenMP):
./build.sh<br>
./feline07.bin 0 1.9 20 7<br>
### detect actual objects and translate into physical properties (redshift, line list) sorted by significance:
python2 mf_new_max_inspect_feline_output.py s2n_v250.fits > catalog.txt<br>
sort -rn -k5 catalog.txt > sorted_catalog.txt<br>
### create comprehensive human readable plots for each detection:
python2 mf_new_plot_and_fit_catalog_allz_wcs_v3.py cube.fits s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028<br>

