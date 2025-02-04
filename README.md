[![status](https://joss.theoj.org/papers/a575acd1ffab0604de7e26eb83fd9bdc/status.svg)](https://joss.theoj.org/papers/a575acd1ffab0604de7e26eb83fd9bdc) [![Docs](http://img.shields.io/badge/Docs-latest-green.svg)](https://feline.readthedocs.io) 


# Feline (Find Emission Lines)
  
## Project Overview

Feline combines a fully parallelized galaxy line template matching with the matched filter approach for individual emission features of LSDcat. 
For the 3D matched filtering, the complete data cube is first median filtered to remove all continuum sources, and then cross-correlated with a template of 
an isolated emission feature in two spatial and one spectral dimension. We assumed a simple Gaussian with a FWHM of $250 \ \rm{km}/\rm{s}$ for the line profile and a PSF 
based on the given seeing in the data cube. The FELINE algorithm then evaluates the likelihood in each spectrum of the cube for emission lines at the positions 
provided by a given redshift and a certain combination of typical emission features. FELINE probes all possible combinations of up to 14 transitions paired 
in 9 groups: $\rm{H}\alpha, \rm{H}\beta, \rm{H}\gamma, \rm{H}\delta$, $\rm{[OII]}$, $\rm{[OIII]}$, $\rm{[NII]}$, $\rm{[SII]}$, and $\rm{[NeIII]}$ for the redshift range of interest $(0.4 < z < 1.4)$. 
This particular selection of lines is motivated by the most prominent emission features expected in the MUSE data within this redshift range. This results in $512 \ (2^9)$ 
different models that are assessed at roughly $8,000$ different redshifts for each of the approx $90,000$ spectra in a single data cube. To ensure that only lines above 
a certain $\rm{S}\rm{N}$ threshold contribute to each model, a penalty value is subtracted for each additional line. The $\rm{S}/\rm{N}$ near strong sky lines are set exactly to that threshold. 
Hence lines that fall onto such a contaminated region will not affect model quality. This is particularly useful for doublet lines that then contribute to a model even when one of 
the lines aligns with a skyline. Furthermore, the $\rm{S}/\rm{N}$ is saturated at a certain threshold to limit the impact of extremely strong lines on the overall budget 
of the tested template. For each spaxel the model with the highest accumulative probability over all contributing lines and its corresponding redshift are determined. 
This approach has the benefit to pick up extremely weak emitters that show multiple emissions lines while avoiding the deluge of false positives when looking 
for single lines below a certain $\rm{S}/\rm{N}$ threshold.

This can be applied to each spatial element independently and was thus fully parallelized. From the resulting spatial map of best model probabilities, 
the peaks were automatically selected via maximum filter and 1D spectra were extracted for each emission line galaxy candidate. The extracted spectra are 
fitted using an emission-line galaxy template, where the redshift and individual line strengths are the only free parameters. This fitting process achieves 
sub-pixel accuracy in the initial redshift estimate while also providing additional diagnostic parameters, such as the $\rm{[OII]}$ doublet ratio, to aid in subsequent 
manual inspection.

## Installation
FELINE requires specific software dependencies for installation and operation. Please follow the instructions below to set up FELINE:

### Prerequisites
To clone the repository with ssh run the following command:
```bash
git clone git@github.com:enthusi/feline.git
```

Ensure the following software is installed on your system:
```bash
python3.x (3.8 or higher)
python3.x-dev
python3.x-venv
clang (recommended due to a significant performance boost compared to gcc) or gcc
SDL2 (Optional: Needed for graphical output during runtime)
```

> [!NOTE] 
> **Mac OS users**: If you use `clang` you only need to install `libomp` e.g. `brew install libomp`. \
> For users which want to use `gcc` only need to adjust the following Makefile lines: \
> `[1] CC = gcc-<version>`  \
> `[2] CFLAGS = -O3 -ffast-math -fopenmp -g -std=c99`
> 
> **Linux users (Debian/Ubuntu)**: If you use `clang` you only need to install `libomp-dev` e.g `apt install libomp-dev`. \
> For users which want to use `gcc` only need to adjust the following Makefile lines: \
> `[1] CC = gcc`  \
> `[2] CFLAGS = -O3 -ffast-math -fopenmp -g -std=c99` 

## Usage Guide
For further information see our [Documentation Website](https://feline.readthedocs.io).

## Contribution Guidelines
We Welcome everyone who wants to Contribute to our Project [Code of Conduct](CODE_OF_CONDUCT.md) and [Contribution Guidelines](CONTRIBUTING.md).


## Acknowledgments
We would like to acknowledge the use of the LSDCat (Line Source Detection and Cataloguing Tool) project for our preprocessing steps. The LSDCat project was developed by Edmund Christian Herenz and has been instrumental in our data analysis.

For more information about LSDCat, please visit the [LSDCat project page](https://bitbucket.org/Knusper2000/lsdcat/src/master/).
