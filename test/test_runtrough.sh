#!/bin/bash

# Test case 1
echo "Running Test Case 1..."
wget http://martinwendt.de/cube.fits -O data/raw/cube.fits --no-verbose  # Download a file
if [ -f data/raw/cube.fits ]; then
  echo "Cube downloaded successfully for Test Case 1"
else
  echo "Failed to download the Cube for Test Case 1"
  echo "Exiting Test on failure (0/5)"
  exit 1 # send error code to github workflow when test failes
fi
echo "____________________________________________________________________"
echo "____________________________________________________________________"

echo "Running Test Case 2..."

#Test Case 2
python src/preprocessing/median-filter-cube.py -t 2 -W 151 --output=data/processed/med_filt.fits data/raw/cube.fits
if [ -f data/processed/med_filt.fits ]; then
  echo "created med_filt.fits successfully for Test Case 2"
else
  echo "Failed to create med_filt.fits for Test Case 2"
  echo "Exiting Test on failure (1/5)"
  exit 1
fi
echo "____________________________________________________________________"
echo "____________________________________________________________________"

echo "Running Test Case 3..."

#Test Case 3
python src/preprocessing/lsd_cc_spatial.py -i data/processed/med_filt.fits -o data/processed/spatial_cc.fits -t 2 -pc 0.7 --lambda0 7050 --gaussian --classic
if [ -f data/processed/spatial_cc.fits ]; then
  echo "created spatial_cc.fits successfully for Test Case 3"
else
  echo "Failed to create spatial_cc.fits for Test Case 3"
  echo "Exiting Test on failure (2/5)"
  exit 1
fi
echo "____________________________________________________________________"
echo "____________________________________________________________________"

echo "Running Test Case 4..."

#Test Case 4
python src/preprocessing/lsd_cc_spectral.py -i data/processed/spatial_cc.fits -F 250 -o data/processed/spectral_cc.fits -t 2 --classic
if [ -f data/processed/spectral_cc.fits ]; then
  echo "created spectral_cc.fits successfully for Test Case 4"
else
  echo "Failed to create spectral_cc.fits for Test Case 4"
  echo "Exiting Test on failure (3/5)"
  exit 1
fi
echo "____________________________________________________________________"
echo "____________________________________________________________________"

#Test Case 5
python src/preprocessing/s2n-cube.py -i data/processed/spectral_cc.fits -o data/processed/s2n_v250.fits -S 1 -N 2 --clobber 
if [ -f data/processed/s2n_v250.fits ]; then
  echo "created s2n_v250.fits successfully for Test Case 5"
else
  echo "Failed to create s2n_v250.fits for Test Case 5"
  echo "Exiting Test on failure (4/5)"
  exit 1
fi
rm data/processed/spectral_cc.fits data/processed/spatial_cc.fits
echo "____________________________________________________________________"
echo "____________________________________________________________________"

echo "####################################################################"
echo "################# Preprocessing succesfull (5/5) ###################"
echo "####################################################################"

exit 0



