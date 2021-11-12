
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR <img src="https://github.com/saezlab/decoupleR/blob/master/inst/figures/logo.svg?raw=1" align="right" width="120" />

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

Many methods allow us to extract biological activities from omics data
using information from prior knowledge resources, reducing the
dimensionality for increased statistical power and better
interpretability. Here, we present decoupleR, a Bioconductor package
containing different statistical methods to extract these signatures
within a unified framework. decoupleR allows the user to flexibly test
any method with any resource. It incorporates methods that take into
account the sign and weight of network interactions. decoupleR can be
used with any omic, as long as its features can be linked to a
biological process based on prior knowledge. For example, in
transcriptomics gene sets regulated by a transcription factor, or in
phospho-proteomics phosphosites that are targeted by a kinase.

<img src="https://github.com/saezlab/decoupleR/blob/master/inst/figures/graphical_abstract.png?raw=1" align="center" width="800">

For more information about how this package has been used with real
data, please check the following links:

-   [decoupleR’s manuscript
    repository](https://github.com/saezlab/decoupleR_manuscript)
-   [Creation of benchmarking
    pipelines](https://github.com/saezlab/decoupleRBench)
-   [Example of Kinase and TF activity
    estimation](https://saezlab.github.io/kinase_tf_mini_tuto/)

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/).

Then install **decoupleR** using from
[Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("decoupleR")

# Check that you have a valid Bioconductor installation
BiocManager::valid()
```

Then install development version from [GitHub](https://github.com/)
with:

``` r
BiocManager::install("saezlab/decoupleR")
```

Or with:

``` r
devtools::install_github("saezlab/decoupleR")
```

## Usage

### Load package and data

We first load the test data included inside **decoupleR**. It consist of
a matrix (`mat`) with logFC coming from transcriptomics, and a
collection of transcription factors that target gene sets with a certain
mode of regulation (either positive or negative).

``` r
library(decoupleR)

inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")

mat <- file.path(inputs_dir, "input-expr_matrix.rds") %>%
    readRDS() %>%
    dplyr::glimpse()
#>  num [1:18490, 1:4] 3.251 0.283 -2.253 0.782 -4.575 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : chr [1:18490] "A1BG" "A1CF" "A2M" "A2ML1" ...
#>   ..$ : chr [1:4] "GSM2753335" "GSM2753336" "GSM2753337" "GSM2753338"

network <- file.path(inputs_dir, "input-dorothea_genesets.rds") %>%
    readRDS() %>%
    dplyr::glimpse()
#> Rows: 151
#> Columns: 5
#> $ tf         <chr> "FOXO4", "FOXO4", "FOXO4", "FOXO4", "FOXO4", "FOXO4", "FOXO…
#> $ confidence <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",…
#> $ target     <chr> "BCL2L11", "BCL6", "CDKN1A", "CDKN1B", "G6PC", "GADD45A", "…
#> $ mor        <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
#> $ likelihood <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,…
```

**Important**: Before running any method in **decoupleR**, we recommend
the user to intersect their prior knowledge with their input matrix
using `intersect_regulons`. This allows to filter out “regulons”
(biological processes) with less than a minimum number of target
features. We recommend to set this value (`minsize`) to at least 5:

``` r
# Remove TFs with less than 5 targets in the input matrix
network <- intersect_regulons(mat, network, tf, target, minsize = 5)
```

### Methods

To check how many methods are currently available in **decoupleR**, run:

``` r
show_methods()
#> # A tibble: 12 × 2
#>    Function      Name                                                           
#>    <chr>         <chr>                                                          
#>  1 run_aucell    AUCell                                                         
#>  2 run_consensus Consensus score between methods                                
#>  3 run_fgsea     Fast Gene Set Enrichment Analysis (FGSEA)                      
#>  4 run_gsva      Gene Set Variation Analysis (GSVA)                             
#>  5 run_mdt       Multivariate Decision Trees (MDT)                              
#>  6 run_mlm       Multivariate Linear Model (MLM)                                
#>  7 run_ora       Over Representation Analysis (ORA)                             
#>  8 run_udt       Univariate Decision Tree (UDT)                                 
#>  9 run_ulm       Univariate Linear Model (ULM)                                  
#> 10 run_viper     Virtual Inference of Protein-activity by Enriched Regulon anal…
#> 11 run_wmean     Weighted Mean (WMEAN)                                          
#> 12 run_wsum      Weighted Sum (WSUM)
```

Function is the function name in **decoupleR** and Name is the full name
of each method. To check how to use individual methods, for example
`mlm`, run `?run_mlm`.

All methods follow the same design pattern and arguments, so moving
between methods should be easy. Here is an example with `mlm`:

``` r
run_mlm(
    mat = mat,
    network = network,
    .source = "tf",
    .target = "target",
    .mor = "mor",
    .likelihood = "likelihood"
)
#> # A tibble: 20 × 5
#>    statistic source condition  score p_value
#>    <chr>     <chr>  <chr>      <dbl>   <dbl>
#>  1 mlm       FOXO4  GSM2753335 2.21  0.0288 
#>  2 mlm       NFIC   GSM2753335 1.12  0.263  
#>  3 mlm       SMAD3  GSM2753335 0.696 0.488  
#>  4 mlm       TFAP2A GSM2753335 1.34  0.183  
#>  5 mlm       RFXAP  GSM2753335 1.63  0.105  
#>  6 mlm       FOXO4  GSM2753336 2.16  0.0325 
#>  7 mlm       NFIC   GSM2753336 1.16  0.250  
#>  8 mlm       SMAD3  GSM2753336 0.882 0.380  
#>  9 mlm       TFAP2A GSM2753336 1.56  0.122  
#> 10 mlm       RFXAP  GSM2753336 2.31  0.0226 
#> 11 mlm       FOXO4  GSM2753337 2.37  0.0195 
#> 12 mlm       NFIC   GSM2753337 0.729 0.467  
#> 13 mlm       SMAD3  GSM2753337 1.11  0.270  
#> 14 mlm       TFAP2A GSM2753337 1.54  0.126  
#> 15 mlm       RFXAP  GSM2753337 2.85  0.00507
#> 16 mlm       FOXO4  GSM2753338 2.10  0.0378 
#> 17 mlm       NFIC   GSM2753338 0.550 0.584  
#> 18 mlm       SMAD3  GSM2753338 0.860 0.391  
#> 19 mlm       TFAP2A GSM2753338 1.28  0.204  
#> 20 mlm       RFXAP  GSM2753338 2.72  0.00742
```

### Decouple wrapper

Moreover, **decoupleR** allows to run multiple methods at the same time
with the function `decouple()`. Statistic functions inside **decoupleR**
always return a tidy tibble that can be easily processed with the tools
provide by the [tidyverse ecosystem](https://www.tidyverse.org/).

``` r
decouple(
    mat = mat,
    network = network,
    .source = "tf",
    .target = "target",
    statistics = c("mlm", "wmean", "ulm", "ora"),
    args = list(
        mlm = list(center=FALSE),
        wmean = list(times=100),
        ulm = list(center=FALSE),
        ora = list(n_up=150, n_bottom=0)
    ),
    consensus_score = TRUE
)
#> # A tibble: 140 × 6
#>    run_id statistic source condition  score p_value
#>     <dbl> <chr>     <chr>  <chr>      <dbl>   <dbl>
#>  1      1 mlm       FOXO4  GSM2753335 2.21   0.0288
#>  2      1 mlm       NFIC   GSM2753335 1.12   0.263 
#>  3      1 mlm       SMAD3  GSM2753335 0.696  0.488 
#>  4      1 mlm       TFAP2A GSM2753335 1.34   0.183 
#>  5      1 mlm       RFXAP  GSM2753335 1.63   0.105 
#>  6      1 mlm       FOXO4  GSM2753336 2.16   0.0325
#>  7      1 mlm       NFIC   GSM2753336 1.16   0.250 
#>  8      1 mlm       SMAD3  GSM2753336 0.882  0.380 
#>  9      1 mlm       TFAP2A GSM2753336 1.56   0.122 
#> 10      1 mlm       RFXAP  GSM2753336 2.31   0.0226
#> # … with 130 more rows
```

It can generate a consensus score between the methods if
`consensus_score = TRUE`.

## Citation

Badia-i-Mompel P., Vélez J., Braunger J., Geiss C., Dimitrov D.,
Müller-Dott S., Taus P., Dugourd A., Holland C.H., Ramirez Flores R.O.
and Saez-Rodriguez J. 2021. decoupleR: Ensemble of computational methods
to infer biological activities from omics data. bioRxiv.
<https://doi.org/10.1101/2021.11.04.467271>

## Contributing to decoupleR

Are you interested in adding a new statistical method or collaborating
in the development of internal tools that allow the extension of the
package? Please check out our [contribution
guide](https://saezlab.github.io/decoupleR/CONTRIBUTING.html).

------------------------------------------------------------------------

Please note that this project is released with a [Contributor Code of
Conduct](https://saezlab.github.io/decoupleR/CODE_OF_CONDUCT). By
participating in this project you agree to abide by its terms.
