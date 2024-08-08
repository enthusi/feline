---
title: 'FELINE: A something...'
tags:
  - Python
  - astronomy
  - ...
authors:
  - name: Martin Wendt
    orcid: 0000-0001-5020-9994
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Marvin Henschel
    orcid: 0009-0000-9462-433X
    affiliation: 2
  - name: Oskar Fjonn Soth
    orcid: 0009-0004-1200-9130
    affiliation: 2

affiliations:
 - name: Institute of Physics and Astronomy, University of Potsdam, 14476 Potsdam, Germany
   index: 1
 - name: Institute of Computer Science and Computational Science, University of Potsdam, 14476 Potsdam, Germany
   index: 2
date: 01 June 2024
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
 1) summary mit allem, aber kurz
 2) science (whatfor and why and how (matched filter), USE CASES
 3) implementation (open MP, C, CUDA, brute force, cache friendly data 
setup), benchmark
 4) references
- github alles in main
- github clean up allgemein
---


# Summary
Feline combines a fully parallelized galaxy line template matching with the matched filter approach for individual emission features of LSDcat [@HerenzE_17a].
For the 3D matched filtering, the complete data cube is first median filtered to remove all continuum sources, and then cross-correlated with a template of an isolated
emission feature in two spatial and one spectral dimension. We assumed a simple Gaussian with a FWHM of 250 km/s for the line profile and a PSF based on the given seeing in the data cube. 
￼
The FELINE algorithm then evaluates the likelihood in each
spectrum of the cube for emission lines at the positions provided by a given redshift and a certain combination of typical emission features. FELINE probes all possible combinations of up to 14 transitions paired in 9 groups: Ha, Hb, Hg, Hd, OII, OIII, NII, SII, and NeIII for the redshift range of interest (0.4 < z < 1.4).
This particular selection of lines is motivated by the most prominent emission features expected in the MUSE data within this redshift range.
This results in 512 (2^9) different models that are assessed at roughly 5,000 different
redshifts for each of the approx 90,000 spectra in a single data cube.
To ensure that only lines above a certain S/N threshold contribute to each
model, a penalty value is subtracted for each additional line. The S/N
near strong sky lines are set exactly to that threshold. Hence lines that fall
onto such a contaminated region will not affect model quality. This is
particularly useful for doublet lines that then contribute to a model even
when one of the lines aligns with a skyline.
Furthermore, the S/N is saturated at a certain threshold to limit the
impact of extremely strong lines on the overall budget of the tested
template. 
For each spaxel the model with the highest accumulative probability over 
all contributing lines and its corresponding redshift are determined.
This approach has the benefit to pick up extremely weak emitters that show
multiple emissions lines while avoiding the deluge of false positives when
looking for single lines below a certain S/N threshold.

This can be applied to each spatial element independently and was thus fully parallelized. From the resulting spatial map of best model probabilities, the peaks were automatically selected via maximum filter and 1D spectra were extracted for each emission line galaxy candidate. Those extracted spectra are fitted with an emission line galaxy template and with the redshift as well as the individual line strengths as only free
parameter to reach sub pixel accuracy in an early redshift estimate as well as deriving further diagnostics for the later manual inspection, such as the OII line ratio.


# Science field

The detection and classification of objects in astrophysical data is a key
task since the very early days of astronomy. Either in image data, or in spectroscopic data, where usually individual emission and absorption signatures are in the focus of interest aside from the underlying continuum fluxes. 
In the past decade the amount of newly observed data has increased dramatically.
In particular the advent of integral field unit spectrographs (IFUs) shifted
the classical single targeted observations by slit spectroscopy to much larger field of views being covered by a single exposure. 
Namely, the VLT/MUSE [@Bacon+10,@2014Msngr.157...13B] 3D spectrograph creates $\sim$ 90,000 medium resolution spectra arranged in a 300 $\times$ 300 spatial grid.
These data cubes have typical sizes of 3-6 GiB per exposure, the sheer amount of data
asks for automated processes to support the scientists.
Simple flux-level peak detection algorithms based on thresholding are prone to either miss alot of potential real objects or as a direct trade off provide a plethora of false-positives.
Fortunately, a couple of assumptions on the expected signals can be made which
enable us to apply algorithms utilizing tis additional information.
One of these approaches is the matched filtering.
The matched filter is the optimal linear filter for maximizing the signal-to-noise ratio (SNR) in the presence of additive stochastic noise, which is a good first approximation of the noise properties we find in our observational data.
The Line Source Detection and Cataloging software [@HerenzE_17a,@herenz2023], a filtering tool developed to recover the single line Ly$\alpha$ emission in surveys with MUSE applies the matched filtering approach.
The resulting signal-to-noise cube reflects the probability of an emission line at the given spatial and spectral position which is significantly boosted by the filtering with a typical line template.
FELINE searches the data cube after the matched-filter was applied for emission features associated with a galaxy.
A physically motivated parameter space of galaxy templates is tested at every spatial position of the cube. Additionally the cosmological redshift $z$ of the galaxies is varied within the observed wavelength range between $0.4 < z < 1.4$.
The detection strength of a galaxy is amplified by the significance for an emission line at the predicted position in the matched filter data cube.
Therefor galaxies with multiple weak emission features can be detected with a 
significance that can strongly exceed that of each individual line.
This approach proves particularly successfull for galaxies that show no continuum flux
in the data and are thus generelly not detected in any imaging data.
We utilized FELINE for our MEGAFLOW survey. The galaxy catalog to which FELINE contributed was explicity used in MF9,10,11.


# Implementation

The FELINE algorithm is based on a brute-force search through the parameter space.
A given data cube of typically 300 x 300 spatial elements, i.e., pixels on the sky
is scanned for all reasonable combinations of a set of lines (including doublets and multiplets) for the whole redshift interval in steps taht correspond to the spectral resolution.
9 different sets of lines were identified, leading to 512 possible combinations.
Each of those is computed for about 5,000 different cosmological redshifts between 0.4 and 1.4 (at this range the prominent OII feature at blabla AA is within the observed wavelength range, which is key for our science case).
This leads to a loop of 300 x 300 x 512 x 5000 = 230,400,000,000 iterations.
Consequently,  C as chosen as the language of implementation for speed.
These numbers also demonstrate the success of this particular approach not to test
full physical models of galaxies (including simulated continuum and temperature broadened emission lines) against the raw observed data, but instead filter the data
with the expected template for a single line and then reducing the individual models to 
a single position at which the likelihood of a line is being probed).
The model has a resolution in wavelength that is limited by the spectral resolution of the data itself and the cosmological redshift only has to be applied to up to 14 discrete values per model rather than the full wavelength range were it to be applied to a full spectral galaxy model.
For each set of parameters (spatial position in the cube, redshift and line composition) the inner FELINE loop returns the value of the highest scoring combination along with its corresponding reshift and line composition. The latter is stored in binary format with a single bit for each line, resulting in a simple 9bit integer number.
This renders the problem ideal for heavy parallelization.
Each of the 300 x 300 spectra in the data cube are considerably small ( < 64KB)
and 512 x 5000 iterations are performed on it with only 3 return values (quality of best match, redshift of best match, line combination of best match).
The outer 300 x 300 iterations are in addition completely independent of each other.

We took advantage of that by full parallelization of that outer loop
in OpenMP with essentially all variables shared due to their independence.
FELINE scales with the number of CPU cores quite well.

Another major improvement in execution time was accomplished by re-arranging the data to maximize the amount of cache hits.
Initially, the raw data is stored as a series images, i.e., 300x300 spatial data points arranged in an array of 4,000 in spectral dimension. The algorithm works on spectral data which caused and interleave of 300x300xsizeof(float) = 360 KB for consecutive data points and the full spectrum exceeding a range of 1 GiB.
As a preprocessing step, the data cube is now re-arranged as a spatial grid of full spectra.

That arrangement plus the fact that individual spectra on which 512 x 5000 iterations are performed are less than 64KB in size
also motivated an implementation of FELINE in CUDA to utilize NVIDA GPUs for parallelization.
Typical full size MUSE data cubes canbe fully loaded into the GPU memory of any modern CUDA capable GPU.
We provide a working implementation that produces identical results to the FELINE C variant.
However, the FELINE C version remains the maintained development version.

Optionally, FELINE plots the three return parameters in real time via SDL surface.
The image dimensions reflect the spatial extent on the sky, the data is color coded following a simple generated color ramp.
Shown are from left to right the quality of the best match, the corresponding redshift of the best match and its template. A fourth panel shows the number of lines that contributed to the most successful model for ease of human readability (it reflects the number of set bits in the best model value).

The data of the three maps is also written to disk.
They contain the full information of the FELINE analysis.

We provide a python framework to further analyse and visualize the FELINE detections.



# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Jon Doe, Jonas Doe, and Doe Jon, and support from Jonas Donas during the genesis of this project.

# References

