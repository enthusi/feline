---
title: 'FELINE: A something...'
tags:
  - Python
  - astronomy
  - ...
authors:
  - name: Martin Wendt
    orcid: 0000-0001-5020-9994
    equal-contrib: true
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Marvin Henschel
    equal-contrib: true
    affiliation: 2
  - name: Oskar Fjonn Soth
    orcid: 0009-0004-1200-9130
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
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
---

# Summary
Feline combines a fully parallelized galaxy line template matching with the matched filter approach for individual emission features of LSDcat `[@HerenzE_17a]`.
For the 3D matched filtering, the complete data cube is first median filtered to remove all continuum sources, and then cross-correlated with a template of an isolated
emission feature in two spatial and one spectral dimension. We assumed a simple Gaussian with a FWHM of 250 km/s for the line profile and a PSF based on the given seeing in the data cube. 
￼
The FELINE algorithm then evaluates the likelihood in each
spectrum of the cube for emission lines at the positions provided by a given redshift and a certain combination of typical emission features. FELINE probes all possible combinations of up to 14 transitions paired in 9 groups: Ha, Hb, Hg, Hd, OII, OIII, NII, SII, and NeIII for the redshift range of interest (0.4 < z < 1.4).
This particular selection of lines is motivated by the most prominent emission features expected in the MUSE data within this redshift range.
This results in 512 (2^9) different models that are assessed at roughly 8,000 different
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


# Statement of need

Here comes the statement of needs


# Mathematics

Here we talk about the math behind our software
Parts of the summary will eventually end up here I guess.
I just wanted to give it a start somehow.

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

