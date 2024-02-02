#!/bin/bash

# Test case 1
echo "Running Test Case 1..."
wget http://martinwendt.de/cube.fits -O data/raw/cube.fits --no-verbose  # Download a file
if [ -f data/raw/cube.fits ]; then
  echo "Cube downloaded successfully for Test Case 1"
else
  echo "Failed to download the Cube for Test Case 1"
  echo "Exiting Test on failure"
  exit 1 # send error code to github workflow when test failes
fi
echo "____________________________________________________________________"
echo "____________________________________________________________________"
echo "Running Test Case 2..."

#Test case 2
python src/preprocessing/median-filter-cube.py data/raw/cube.fits --signalHDU=1 --varHDU=2 --num_cpu=2 --width=151 --output=data/processed/med_filt.fits

if [-f data/processed/med_filt.fits]; then
  echo "created med_filt.cubes successfully for Test Case 2"
else
  echo "Failed to create med_filt.cubes for Test Case 2"
  echo "Exiting Test on failure"
  exit 1
fi

exit 0



