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
  - name: Oskar Fjonn Soth
    orcid: 0009-0004-1200-9130
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Marvin Henschel
    equal-contrib: true
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

Here the software summary belongs

# Statement of need

Here comes the statement of needs

# Mathematics

Here we talk about the math behind our software

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
