
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/saezlab/decoupleR.svg?token=PagY1pyvMyyL3AJHRy5V&branch=master)](https://travis-ci.com/saezlab/decoupleR)
[![Lifecycle:experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/saezlab/decoupleR/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/decoupleR?branch=master)
<!-- badges: end -->

> a community effort by saezlab members

## Overview

**Under development** - This R package allows to combine a variety of
gene sets with a variety of statistics for functional genomics analyses.

## Before you start

### Dependencies

To ensure that we all use the same package versions we make use of
[renv](https://rstudio.github.io/renv/articles/renv.html), a package for
dependency management.

``` r
# install all required packages using the renv package
renv::restore()
```

## How to add a new statistic

  - Open an [issue](https://github.com/saezlab/decoupleR/issues) stating
    which statistic you would like to add. Assign it to yourself and
    label it with the badge `enhancement`. Please check before if not
    someone else is already working on implementing this statistic.

  - Don’t work on the `master branch`. Either create a new `branch` with
    a meaningful name or make a `fork`. When the function is implemented
    and tested (\!) we will use `pull requests` to integrate the new
    feature in the `master branch`.

  - Assuming you would like to implement `GSEA`. You will need to define
    the following set of functions:
    
      - All following functions will be written in the script
        `R/gsea.R`.
      - `run_gsea(emat, genesets, list=options(), gs_resource, tidy)`.
        Please check the function
        [run\_viper](https://github.com/saezlab/decoupleR/blob/master/R/viper.R#L30)
        to understand the arguments.
      - `make_gsea_genesets(genesets)`. First define a standardized
        input for gene sets in *tibble/dataframe* format for your
        statistic. Then this function should convert this table in the
        required input for the underlying function (in this case list of
        gene sets for the function `fgsea()`). Check an example
        [here](https://github.com/saezlab/decoupleR/blob/master/R/viper.R#L68).
        In the roxygen comments of the function you define which columns
        must be available in the standardized format.
      - In parallel you need to define helper functions that convert
        every available gene set (e.g. progeny and DoRothEA gene sets)
        to your defined standardized input (e.g. `progeny2gsea()`).
        Check an example
        [here](https://github.com/saezlab/decoupleR/blob/master/R/viper.R#L90)
      - **All functions must be documented following the roxygen2
        standard.**

  - We will use unit tests for the `run_x()` functions to ensure that
    our functions work properly.
    [Here](https://github.com/saezlab/decoupleR/blob/master/tests/testthat/test-viper.R)
    are examples for the `run_viper` function.

  - *Optional*: For a consistent coding style and efficient
    implementation we will mainly use `tidyverse`.

## How to add a new gene set resource

  - Open an [issue](https://github.com/saezlab/decoupleR/issues) stating
    which gene set resource you would like to add. Assign it to yourself
    and label it with the badge `enhancement`. Please check before if
    not someone else is already working on implementing this gene set
    resource.

  - Don’t work on the `master branch`. Either create a new `branch` with
    a meaningful name or make a `fork`. When the function is implemented
    and tested (\!) we will use `pull requests` to integrate the new
    feature in the `master branch`.

  - Define helper functions for each available statistics that convert
    your gene sets to the standardized format of the respective
    statistic (e.g. `your_genesets2viper()`).

  - Deposit a representative selection of gene sets in the directory
    `inst/testdata` and implement unit tests for all available
    statistics coupled with your new gene sets.

## Installation

``` r
# install the development version from GitHub
# install.packages("remotes")
remotes::install_github("saezlab/decoupleR")
```
