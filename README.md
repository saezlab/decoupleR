
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR

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

Computational methods allow the extraction of mechanistic signatures
from omics data based on prior knowledge resources, reducing the
dimensionality of the data for increased statistical power and better
interpretability. Here, we present decoupleR, a Bioconductor package
containing different statistical methods to extract these signatures
within a unified framework. decoupleR allows the user to flexibly test
any method with any resource. It incorporates methods that take into
account the sign and weight of network interactions.

For more information about how this package has been used with real
data, please check the following links:

-   [decoupleR’s manuscript
    repository](https://github.com/saezlab/decoupleR_manuscript)
-   [Creation of benchmarking
    pipelines](https://github.com/saezlab/decoupleRBench).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/).

Then install `decoupleR` using from
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

## Usage

### Load packaga and data

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

### Decouple wrapper

`decouple()` allows access to all **decoupleR** available statistics in
one place. Statistic functions inside **decoupleR** always return a tidy
tibble that can be easily processed with the tools provide by the
[tidyverse ecosystem](https://www.tidyverse.org/).

``` r
decouple(
    mat = mat,
    network = network,
    .source = "tf",
    .target = "target",
    statistics = c("gsva", "wmean", "wsum", "ulm", "ora"),
    args = list(
        gsva = list(verbose = FALSE),
        wmean = list(times=100),
        wsum = list(times=100),
        ulm = list(center=FALSE),
        ora = list()
    )
)
#> # A tibble: 200 × 6
#>    run_id statistic source condition    score p_value
#>    <chr>  <chr>     <chr>  <chr>        <dbl>   <dbl>
#>  1 1      gsva      FOXO4  GSM2753335 -0.380       NA
#>  2 1      gsva      FOXO4  GSM2753336 -0.300       NA
#>  3 1      gsva      FOXO4  GSM2753337  0.239       NA
#>  4 1      gsva      FOXO4  GSM2753338  0.0907      NA
#>  5 1      gsva      NFIC   GSM2753335 -0.0845      NA
#>  6 1      gsva      NFIC   GSM2753336  0.0778      NA
#>  7 1      gsva      NFIC   GSM2753337 -0.260       NA
#>  8 1      gsva      NFIC   GSM2753338  0.281       NA
#>  9 1      gsva      RFXAP  GSM2753335 -0.810       NA
#> 10 1      gsva      RFXAP  GSM2753336 -0.472       NA
#> # … with 190 more rows
```

### Individual parts

In turn, we recognize that the use of individual statistics may be of
interest. Therefore, these are also exported and ready for use. All
statistics follow the same design pattern and arguments, so moving
between statistics could be very comfortable.

``` r
# wmean call is equivalent to the one made by decouple() above.
run_wmean(
    mat = mat,
    network = network,
    .source = "tf",
    .target = "target",
    .likelihood = "likelihood"
)
#> # A tibble: 60 × 5
#>    statistic  source condition  score p_value
#>    <chr>      <chr>  <chr>      <dbl>   <dbl>
#>  1 corr_wmean FOXO4  GSM2753335  5.94    0   
#>  2 corr_wmean FOXO4  GSM2753336  5.97    0   
#>  3 corr_wmean FOXO4  GSM2753337  6.36    0   
#>  4 corr_wmean FOXO4  GSM2753338  6.40    0   
#>  5 corr_wmean NFIC   GSM2753335  1.58    0.23
#>  6 corr_wmean NFIC   GSM2753336  1.84    0.2 
#>  7 corr_wmean NFIC   GSM2753337  1.23    0.27
#>  8 corr_wmean NFIC   GSM2753338  1.27    0.28
#>  9 corr_wmean RFXAP  GSM2753335  3.61    0.11
#> 10 corr_wmean RFXAP  GSM2753336  6.51    0.04
#> # … with 50 more rows
```

<!-- ## Citation -->
<!-- Below is the citation output from using `citation('decoupleR')` in R. Please -->
<!-- run this yourself to check for any updates on how to cite __decoupleR__. -->
<!-- ```{r 'citation', eval = requireNamespace('decoupleR')} -->
<!-- print(citation("decoupleR"), bibtex = TRUE) -->
<!-- ``` -->
<!-- Please note that the `decoupleR` was only made possible thanks to many other R -->
<!-- and bioinformatics software authors, which are cited either in the vignettes -->
<!-- and/or the paper(s) describing this package. -->

## Contributing to decoupleR

Are you interested in adding a new statistical method or collaborating
in the development of internal tools that allow the extension of the
package? Please check out our [contribution
guide](https://saezlab.github.io/decoupleR/CONTRIBUTING.html).

------------------------------------------------------------------------

Please note that this project is released with a [Contributor Code of
Conduct](https://saezlab.github.io/decoupleR/CODE_OF_CONDUCT). By
participating in this project you agree to abide by its terms.
