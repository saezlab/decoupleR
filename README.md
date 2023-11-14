
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR <img src="inst/figures/logo.svg" align="right" width="120" />

<!-- badges: start -->

[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![BioC
status](http://www.bioconductor.org/shields/build/release/bioc/decoupleR.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/decoupleR)
[![BioC dev
status](http://www.bioconductor.org/shields/build/devel/bioc/decoupleR.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/decoupleR)
[![R build
status](https://github.com/saezlab/decoupleR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/saezlab/decoupleR/actions)
[![Codecov test
coverage](https://codecov.io/gh/saezlab/decoupleR/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/decoupleR?branch=master)
[![GitHub
issues](https://img.shields.io/github/issues/saezlab/decoupleR)](https://github.com/saezlab/decoupleR/issues)
<!-- badges: end -->

## Overview

There are many methods that allow us to extract biological activities
from omics data. `decoupleR` is a Bioconductor package containing
different statistical methods to extract biological signatures from
prior knowledge within a unified framework. Additionally, it
incorporates methods that take into account the sign and weight of
network interactions. `decoupleR` can be used with any omic, as long as
its features can be linked to a biological process based on prior
knowledge. For example, in transcriptomics gene sets regulated by a
transcription factor, or in phospho-proteomics phosphosites that are
targeted by a kinase. This is the R version, for its faster and memory
efficient Python implementation go
[here](https://decoupler-py.readthedocs.io/en/latest/).

<p align="center" width="100%">
<img src="https://github.com/saezlab/decoupleR/blob/master/inst/figures/graphical_abstract.png?raw=1" align="center" width="45%">
</p>

For more information about how this package has been used with real
data, please check the following links:

- [decoupleR’s general
  usage](https://saezlab.github.io/decoupleR/articles/decoupleR.html)
- [Pathway activity inference in bulk
  RNA-seq](https://saezlab.github.io/decoupleR/articles/pw_bk.html)
- [Pathway activity inference from
  scRNA-seq](https://saezlab.github.io/decoupleR/articles/pw_sc.html)
- [Transcription factor activity inference in bulk
  RNA-seq](https://saezlab.github.io/decoupleR/articles/tf_bk.html)
- [Transcription factor activity inference from
  scRNA-seq](https://saezlab.github.io/decoupleR/articles/tf_sc.html)
- [Example of Kinase and TF activity
  estimation](https://saezlab.github.io/kinase_tf_mini_tuto/)
- [decoupleR’s manuscript
  repository](https://github.com/saezlab/decoupleR_manuscript)
- [Python
  implementation](https://decoupler-py.readthedocs.io/en/latest/)

# Installation

`decoupleR` is an R package distributed as part of the Bioconductor
project. To install the package, start R and enter:

``` r
install.packages("BiocManager")
BiocManager::install("saezlab/decoupleR")
```

## License

Footprint methods inside `decoupleR` can be used for academic or
commercial purposes, except `viper` which holds a non-commercial
license.

The data redistributed by `OmniPath` does not have a license, each
original resource carries their own. [Here](https://omnipathdb.org/info)
one can find the license information of all the resources in `OmniPath`.

## Citation

Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov
D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores
R.O. and Saez-Rodriguez J. 2022. decoupleR: ensemble of computational
methods to infer biological activities from omics data. Bioinformatics
Advances. <https://doi.org/10.1093/bioadv/vbac016>
