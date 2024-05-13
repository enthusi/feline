#!/bin/bash

# Delete all files except .gitkeep in data/raw
find data/raw ! -name '.gitkeep' -type f -delete

# Delete all files except .gitkeep in data/processed
find data/processed ! -name '.gitkeep' -type f -delete

# Remove all files with .bmp and .raw extension from the current directory
rm -f *.bmp *.raw

# Delete all files with .txt and .fits extension in src/postprocessing
find src/postprocessing -type f \( -name '*.txt' -o -name '*.fits' -o -name '*.png' -o -name '*.log' -o ! -name "*.*" \) -delete



