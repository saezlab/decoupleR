
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/saezlab/decoupleR.svg?token=PagY1pyvMyyL3AJHRy5V&branch=master)](https://travis-ci.com/saezlab/decoupleR)
[![Lifecycle:experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/saezlab/decoupleR/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/decoupleR?branch=master)
<!-- badges: end -->

> a community effort by [saezlab](http://saezlab.org) members

## Overview

**Under development** - decoupleR aims to combine various gene sets
resources with a variety of statistics for functional genomics analyses.

## Installation

``` r
# install the development version from GitHub
# install.packages("devtools")
devtools::install_github("saezlab/decoupleR")
devtools::install_github("saezlab/decoupleR@devel-jesus")
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
#> $ tf         <chr> "FOXO4", "FOXO4", "FOXO4", "FOXO4", "FOXO4", "FOXO4", "FOX…
#> $ confidence <chr> "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"…
#> $ target     <chr> "BCL2L11", "BCL6", "CDKN1A", "CDKN1B", "G6PC", "GADD45A", …
#> $ mor        <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
#> $ likelihood <dbl> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1…
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
  statistics = c("gsva", "mean", "pscira", "scira", "viper"),
  args = list(
    gsva = list(verbose = FALSE),
    mean = list(.mor = "mor", .likelihood = "likelihood"),
    pscira = list(.mor = "mor"),
    scira = list(.mor = "mor"),
    viper = list(.mor = "mor", .likelihood = "likelihood", verbose = FALSE)
  )
)
#> # A tibble: 112 x 7
#>    run_id statistic tf    condition    score statistic_time p_value
#>    <chr>  <chr>     <chr> <chr>        <dbl> <drtn>           <dbl>
#>  1 1      gsva      FOXO4 GSM2753335 -0.380  6.696238 secs       NA
#>  2 1      gsva      FOXO4 GSM2753336 -0.300  6.696238 secs       NA
#>  3 1      gsva      FOXO4 GSM2753337  0.239  6.696238 secs       NA
#>  4 1      gsva      FOXO4 GSM2753338  0.0907 6.696238 secs       NA
#>  5 1      gsva      NFIC  GSM2753335 -0.0845 6.696238 secs       NA
#>  6 1      gsva      NFIC  GSM2753336  0.0778 6.696238 secs       NA
#>  7 1      gsva      NFIC  GSM2753337 -0.260  6.696238 secs       NA
#>  8 1      gsva      NFIC  GSM2753338  0.281  6.696238 secs       NA
#>  9 1      gsva      RFXAP GSM2753335 -0.810  6.696238 secs       NA
#> 10 1      gsva      RFXAP GSM2753336 -0.472  6.696238 secs       NA
#> # … with 102 more rows
```

### Individual parts

In turn, we recognize that the use of individual statistics may be of
interest. Therefore, these are also exported and ready for use. All
statistics follow the same design pattern and arguments, so moving
between statistics could be very comfortable.

``` r
# viper call is equivalent to the one made by decouple() above.
run_viper(
  mat = mat,
  network = network,
  .source = "tf",
  .target = "target",
  .likelihood = "likelihood",
  verbose = FALSE
)
#> # A tibble: 12 x 5
#>    statistic tf     condition    score statistic_time
#>    <chr>     <chr>  <chr>        <dbl> <drtn>        
#>  1 viper     NFIC   GSM2753335  0.0696 0.0346725 secs
#>  2 viper     NFIC   GSM2753336 -0.0265 0.0346725 secs
#>  3 viper     NFIC   GSM2753337 -0.516  0.0346725 secs
#>  4 viper     NFIC   GSM2753338 -0.543  0.0346725 secs
#>  5 viper     SMAD3  GSM2753335  0.176  0.0346725 secs
#>  6 viper     SMAD3  GSM2753336  0.0426 0.0346725 secs
#>  7 viper     SMAD3  GSM2753337  0.219  0.0346725 secs
#>  8 viper     SMAD3  GSM2753338  0.142  0.0346725 secs
#>  9 viper     TFAP2A GSM2753335  0.722  0.0346725 secs
#> 10 viper     TFAP2A GSM2753336  0.582  0.0346725 secs
#> 11 viper     TFAP2A GSM2753337  0.462  0.0346725 secs
#> 12 viper     TFAP2A GSM2753338  0.330  0.0346725 secs
```

## Contributing to decoupleR

Are you interested in adding a new statistical method or collaborating
in the development of internal tools that allow the extension of the
package? Please check out our [contribution
guide](https://saezlab.github.io/decoupleR/CONTRIBUTING.html).

-----

Please note that this project is released with a [Contributor Code of
Conduct](https://saezlab.github.io/decoupleR/CODE_OF_CONDUCT). By
participating in this project you agree to abide by its terms.
