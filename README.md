# Feline (Find Emission Lines)
## Table of Contents
- [Project Overview](#project-overview)
- [Installation Instructions](#installation-instructions)
- [Usage Guide](#usage-guide)
- [Features](#features)
- [Documentation](#documentation)
- [Contribution Guidelines](#contribution-guidelines)
- [License](#license)
- [Changelog](#changelog)
- [Contact Information](#contact-information)
- [Acknowledgments](#acknowledgments)
  
## Project Overview
## Project Overview

Feline combines a fully parallelized galaxy line template matching with the matched filter approach for individual emission features of LSDcat. For the 3D matched filtering, the complete data cube is first median filtered to remove all continuum sources, and then cross-correlated with a template of an isolated emission feature in two spatial and one spectral dimension. We assumed a simple Gaussian with a FWHM of 250 km/s for the line profile and a PSF based on the given seeing in the data cube. ￼ The FELINE algorithm then evaluates the likelihood in each spectrum of the cube for emission lines at the positions provided by a given redshift and a certain combination of typical emission features. FELINE probes all possible combinations of up to 14 transitions paired in 9 groups: Ha, Hb, Hg, Hd, OII, OIII, NII, SII, and NeIII for the redshift range of interest (0.4 < z < 1.4). This particular selection of lines is motivated by the most prominent emission features expected in the MUSE data within this redshift range. This results in 512 (2^9) different models that are assessed at roughly 8,000 different redshifts for each of the approx 90,000 spectra in a single data cube. To ensure that only lines above a certain S/N threshold contribute to each model, a penalty value is subtracted for each additional line. The S/N near strong sky lines are set exactly to that threshold. Hence lines that fall onto such a contaminated region will not affect model quality. This is particularly useful for doublet lines that then contribute to a model even when one of the lines aligns with a skyline. Furthermore, the S/N is saturated at a certain threshold to limit the impact of extremely strong lines on the overall budget of the tested template. For each spaxel the model with the highest accumulative probability over all contributing lines and its corresponding redshift are determined. This approach has the benefit to pick up extremely weak emitters that show multiple emissions lines while avoiding the deluge of false positives when looking for single lines below a certain S/N threshold.

This can be applied to each spatial element independently and was thus fully parallelized. From the resulting spatial map of best model probabilities, the peaks were automatically selected via maximum filter and 1D spectra were extracted for each emission line galaxy candidate. Those extracted spectra are fitted with an emission line galaxy template and with the redshift as well as the individual line strengths as only free parameter to reach sub pixel accuracy in an early redshift estimate as well as deriving further diagnostics for the later manual inspection, such as the OII line ratio.

## Installation Instructions
The following software is required to work with the software
```bash
python3.x (3.8 or higher)
python3.x-dev
python3.x-venv
gcc or clang
SDL2 (Optional: Only necessary when you want to have graphical output during the run)
```
## Usage Guide
You can either excecute the program manually or with our prewritten workflow inside the make file
### Automatic Usage
Edit the CUBENAME and CUBELINK inside the Makefile

![image](https://github.com/user-attachments/assets/d6f2383a-e3e2-4a55-910e-9401bfe3cbea)

If you have your CUBEFILE stored locally copy the cube in the root folder of the Project and just edit the CUBENAME.

After editing the Makefile you can run the workflow with
```bash
make run
```
When the Workflow finished and all PDF Files are successfully merged you can find you results at
```bash
data/pdf_files/result_*.pdf
```
You can clean up all temporary files with
```bash
make clean
```
### Manual Usage
Copy the Cubefile into the data/raw/ directory in the project

Create a virtual environment
```bash
python3.8 -m venv venv
```
Activate the environment
```bash
source venv/bin/activate
```
Install all required python packages
```bash
pip install -r requirements.txt
```
Now we can start with the preprocessing

First we apply a median filter to 'flat' out the data
```bash
python src/preprocessing/median-filter-cube.py data/raw/<CUBENAME>.fits --signalHDU=1 --varHDU=2 --num_cpu=<num_cores> --width=151 --output=data/processed/med_filt.fits
```
now we filter the data cube with a 'typical line' template in spatial dimension
```bash
python src/preprocessing/lsd_cc_spatial.py --input=data/processed/med_filt.fits --SHDU=1 --NHDU=2 --threads=<num_cores> --gaussian --lambda0=7050 -pc 0.7 --classic --output=data/processed/spatial_cc.fits --overwrite
```
now we filter the data cube with a 'typical line' template in spectral dimension
```bash
python src/preprocessing/lsd_cc_spectral.py --input=data/processed/spatial_cc.fits --threads=<num_cores> --FWHM=250 --SHDU=1 --NHDU=2 --classic --output=data/processed/spectral_cc.fits --overwrite
```
And finally we construct a signal-to-noise cube
```bash
python src/preprocessing/s2n-cube.py --input=data/processed/spectral_cc.fits --output=data/processed/s2n_v250.fits --clobber --NHDU=2 --SHDU=1
```
Before we start with the Main Program we will transpose the Cube for better cache access
```
cd src/postprocessing
python combination.py <cubename>.fits s2n_v250.fits
```
Now we can compile and run the main program
```bash
cd ../../
make
./ feline.bin ZLOW ZHIGH MAX_MATCH IGNORE_BELOW
```
Now we start the postprocessing
```bash
cd src/postprocessing
```
we start by running a script to detect all objects inside the signale-to-noise cube and safe the results in a sorted_catalog.txt file
```bash
python detect_objects.py s2n_v250.fits
```
next we will create the final plots for our found objects and merge the pdf files into one
```bash
python create_final_plots.py <cubename>.fits s2n_v250.fits sorted_catalog.txt med_filt.fits J0014m0028
python create_pdf.py
```
You can find your results in
```bash
data/pdf_files/result_*.pdf
```
after running the workflow you can clean up all temporary files with
```bash
make clean
```
## Features
## Documentation
## Contribution Guidelines
We Welcome everyone who wants to Contribute to our Project you can find our Code of Conduct and Contribution Guidelines here:

[![Contributor Covenant](https://img.shields.io/badge/Contributor_Convenant-2.1-blue)](CONDUCT.md)
[![CONTRIBUTION](https://img.shields.io/badge/CONTRIBUTING.md-V1.0-blue)](CONTRIBUTING.md)
## License
## Changelog
## Contact Information
## Acknowledgments
We would like to acknowledge the use of the LSDCat (Line Source Detection and Cataloguing Tool) project for our preprocessing steps. The LSDCat project was developed by Edmund Christian Herenz and has been instrumental in our data analysis.

For more information about LSDCat, please visit the [LSDCat project page](https://bitbucket.org/Knusper2000/lsdcat/src/master/).
