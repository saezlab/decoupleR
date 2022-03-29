
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
targeted by a kinase.

<p align="center" width="100%">
<img src="https://github.com/saezlab/decoupleR/blob/master/inst/figures/graphical_abstract.png?raw=1" align="center" width="45%">
</p>

For more information about how this package has been used with real
data, please check the following links:

-   [decoupleR’s
    vignette](https://saezlab.github.io/decoupleR/articles/decoupleR.html)
-   [Python
    implementation](https://decoupler-py.readthedocs.io/en/latest/)
-   [decoupleR’s manuscript
    repository](https://github.com/saezlab/decoupleR_manuscript)
-   [Creation of benchmarking
    pipelines](https://github.com/saezlab/decoupleRBench)
-   [Example of Kinase and TF activity
    estimation](https://saezlab.github.io/kinase_tf_mini_tuto/)

# Installation

`decoupleR` is an R package distributed as part of the Bioconductor
project. To install the package, start R and enter:

``` r
install.packages("BiocManager")
BiocManager::install("decoupleR")
```

Alternatively, you can instead install the latest development version
from [GitHub](https://github.com/) with:

``` r
BiocManager::install("saezlab/decoupleR")
```

## Citation

Badia-i-Mompel P., Vélez Santiago J., Braunger J., Geiss C., Dimitrov
D., Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores
R.O. and Saez-Rodriguez J. 2022. decoupleR: ensemble of computational
methods to infer biological activities from omics data. Bioinformatics
Advances. <https://doi.org/10.1093/bioadv/vbac016>
