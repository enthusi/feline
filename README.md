# feline
in parent dir:<br>
wget http://martinwendt.de/cube.fits

in preprocess directory:<br>
python2 median-filter-cube.py ../cube.fits --signalHDU=1 --varHDU=2 --num_cpu=6 --width=151 --output=med_filt.fits<br>
python2 lsd_cc_spatial_mask.py --input=med_filt.fits --SHDU=0 --NHDU=1 --threads=4 --gaussian --lambda0=7050 -p0=0.7 --output=spatial_cc.fits<br>
python2 lsd_cc_spectral.py --input=spatial_cc.fits --threads=2 --FWHM=250 --SHDU=0 --NHDU=1 --output=spectral_cc.fits<br>
python2 s2n-cube.py --input=spectral_cc.fits --output=s2n_v250.fits --clobber<br>

copy s2n_v250.fits and med_filt.fits into parent dir<br>

in parent dir:<br>
python2 createnan.py cube.fits<br>
python2 preprocess_s2n_cube.py s2n_v250.fits none<br>
./build.sh<br>
./feline07.bin 0 1.9 20 7<br>
python2 mf_new_max_inspect_feline_output.py s2n_v250.fits > catalog.txt<br>
sort -rn -k5 catalog.txt > sorted_catalog.txt<br>
python2 mf_new_plot_and_fit_catalog_allz_wcs_v3.py cube.fits s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028<br>

