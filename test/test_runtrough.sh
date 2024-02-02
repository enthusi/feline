#!/bin/bash

# Test case 1
echo "Running Test Case 1..."
wget https://martinwendt.de/cube.fits  # Download a file
if [ -f cube.fits ]; then
  echo "Cube downloaded successfully for Test Case 1"
else
  echo "Failed to download the Cube for Test Case 1"
  exit 1 # send error code to github workflow when test failes
fi


