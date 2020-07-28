
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR

<!-- badges: start -->

[![Build
Status](https://travis-ci.com/saezlab/decoupleR.svg?token=PagY1pyvMyyL3AJHRy5V&branch=master)](https://travis-ci.com/saezlab/decoupleR)
[![Lifecycle:experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Codecov test
coverage](https://codecov.io/gh/saezlab/decoupleR/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/decoupleR?branch=master)
<!-- badges: end -->

> a community effort by [saezlab](www.saezlab.org) members

## Overview

**Under development** - decoupleR aims to combine various gene sets
resources with a variety of statistics for functional genomics analyses.

## How to contribute?

### How to add a new statistic

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
        every available gene set (e.g. progeny and dorothea gene sets)
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

### How to add a new gene set resource

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

### How to integrate your changes in the master branch

Before you push your changes to your development branch please first
make sure that `devtools::check()` runs without any problems on your
local machine (warning(s) might be acceptable, but will be decided on a
case-by-case basis). After you pushed then to your development branch
[travis](https://travis-ci.com/github/saezlab/decoupleR) will also build
(`R CMD build`) and check (`R CMD check`) the package. **Your changes
will only be integrated in the master branch when the package passes all
checks on travis.**

If all requirements are fulfilled (passed `R CMD check` on your local
machine and travis, implemented unit tests, proper documentation), add a
brief and concise summary of the implemented features in the
[NEWS.md](https://github.com/saezlab/decoupleR/blob/master/NEWS.md) file
and bump the version in the
[DESCRIPTION](https://github.com/saezlab/decoupleR/blob/master/DESCRIPTION)
file. In addition list you as contributor (ctr) / author (aut) also in
the
[DESCRIPTION](https://github.com/saezlab/decoupleR/blob/master/DESCRIPTION)
file.

After you open a pull request please refer in the commit message to the
related issue, e.g. `Closes #5`. This message automatically closes issue
5, after successful integration. Also request review from someone else
to double check your code. **Do not accept your own pull request without
someone else having checked your code**.

## Installation

``` r
# install the development version from GitHub
# install.packages("remotes")
remotes::install_github("saezlab/decoupleR")
```
