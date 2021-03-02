
<!-- README.md is generated from README.Rmd. Please edit that file -->

# decoupleR

<!-- badges: start -->

[![Lifecycle:experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![R build
status](https://github.com/saezlab/decoupleR/workflows/R-CMD-check-bioc/badge.svg)](https://github.com/saezlab/decoupleR/actions)
[![Codecov test
coverage](https://codecov.io/gh/saezlab/decoupleR/branch/master/graph/badge.svg)](https://codecov.io/gh/saezlab/decoupleR?branch=master)
[![GitHub
issues](https://img.shields.io/github/issues/saezlab/decoupleR)](https://github.com/saezlab/decoupleR/issues)
<!-- badges: end -->

<!-- > A community effort by [saezlab](http://saezlab.org) members. -->

## Overview

Transcriptome profiling followed by differential gene expression
analysis often leads to lists of genes that are hard to analyze and
interpret. Downstream analysis tools can be used to summarize
deregulation events into a smaller set of biologically interpretable
features. In particular, methods that estimate the activity of
transcription factors (TFs) from gene expression are commonly used. It
has been shown that the transcriptional targets of a TF yield a much
more robust estimation of the TF activity than observing the expression
of the TF itself. Consequently, for the estimation of transcription
factor activities, a network of transcriptional regulation is required
in combination with a statistical algorithm that summarizes the
expression of the target genes into a single activity score. Over the
years, many different regulatory networks and statistical algorithms
have been developed, mostly in a fixed combination of one network and
one algorithm. To systematically evaluate both networks and algorithms,
we developed decoupleR , an R package that allows users to apply
efficiently any combination provided.

As an initial set of regulatory networks, we integrated the following
resources:

-   [DoRothEA](https://github.com/saezlab/dorothea)
-   [CHEA3](https://amp.pharm.mssm.edu/ChEA3)
-   [RegNetwork](http://www.regnetworkweb.org/)

And the following TF activity inference methods:

### Estimation based on enrichment of sets on rankings

-   [viper](http://bioconductor.org/packages/release/bioc/html/viper.html)
-   GSEA as implemented in
    [fgsea](https://www.bioconductor.org/packages/release/bioc/html/fgsea.html)
-   [GSVA](https://www.bioconductor.org/packages/release/bioc/html/GSVA.html)

### Estimation from linear models

-   SCIRA: Linear models to estimate TF activities from gene expression
    as defined
    [here](https://www.nature.com/articles/s41525-020-00151-y?elqTrackId=d7efb03cf5174fe2ba84e1c34d602b13)

-   pscira: Linear combination of gene expression based on mode of
    regulation followed by a comparison to a random null model.

### Estimation from statistics

-   mean: Weighted mean that allows the use of directions and
    contribution weights.

-   normalized\_mean: Similar as above, however the final score is
    corrected based on a null (random) model

We benchmarked an initial set of 84 combinations, comprising 7 methods
and 13 networks.

We evaluated the precision of different combinations in recovering
perturbed TFs from different collections of gene expression datasets.
Additionally, we tested the effects of combining multiple sources and
estimations. We set up the package in a modular way which makes it easy
and intuitive to extend it with further statistics or networks. We
invite the community to participate by implementing their own statistics
or integrating their gene regulatory network. With the decoupleR
package, we lay the foundation for a crowdsourced systematic assessment
of transcription factor activity estimation from transcriptomics data.

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/).

<!-- Then install `decoupleR` using from [Bioconductor](http://bioconductor.org/) the following code: -->
<!-- ```{r bioconductor_install, eval = FALSE} -->
<!-- if (!requireNamespace("BiocManager", quietly = TRUE)) { -->
<!--     install.packages("BiocManager") -->
<!-- } -->
<!-- BiocManager::install("decoupleR") -->
<!-- ``` -->

Then install development version from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("saezlab/decoupleR")
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
    statistics = c("gsva", "mean", "pscira", "scira", "viper", "ora"),
    args = list(
        gsva = list(verbose = FALSE),
        mean = list(.mor = "mor", .likelihood = "likelihood"),
        pscira = list(.mor = "mor"),
        scira = list(.mor = "mor"),
        viper = list(.mor = "mor", .likelihood = "likelihood", verbose = FALSE),
        ora = list()
    )
)
#> Warning: multiple methods tables found for 'rowRanges'
#> # A tibble: 140 x 11
#>    run_id statistic tf    condition   score p_value estimate conf.low conf.high
#>    <chr>  <chr>     <chr> <chr>       <dbl>   <dbl>    <dbl>    <dbl>     <dbl>
#>  1 1      gsva      FOXO4 GSM27533… -0.380       NA       NA       NA        NA
#>  2 1      gsva      FOXO4 GSM27533… -0.300       NA       NA       NA        NA
#>  3 1      gsva      FOXO4 GSM27533…  0.239       NA       NA       NA        NA
#>  4 1      gsva      FOXO4 GSM27533…  0.0907      NA       NA       NA        NA
#>  5 1      gsva      NFIC  GSM27533… -0.0845      NA       NA       NA        NA
#>  6 1      gsva      NFIC  GSM27533…  0.0778      NA       NA       NA        NA
#>  7 1      gsva      NFIC  GSM27533… -0.260       NA       NA       NA        NA
#>  8 1      gsva      NFIC  GSM27533…  0.281       NA       NA       NA        NA
#>  9 1      gsva      RFXAP GSM27533… -0.810       NA       NA       NA        NA
#> 10 1      gsva      RFXAP GSM27533… -0.472       NA       NA       NA        NA
#> # … with 130 more rows, and 2 more variables: method <chr>, alternative <chr>
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
#> # A tibble: 20 x 4
#>    statistic tf     condition    score
#>    <chr>     <chr>  <chr>        <dbl>
#>  1 viper     FOXO4  GSM2753335  1.34  
#>  2 viper     FOXO4  GSM2753336  1.18  
#>  3 viper     FOXO4  GSM2753337  1.44  
#>  4 viper     FOXO4  GSM2753338  1.21  
#>  5 viper     NFIC   GSM2753335  0.0696
#>  6 viper     NFIC   GSM2753336 -0.0265
#>  7 viper     NFIC   GSM2753337 -0.516 
#>  8 viper     NFIC   GSM2753338 -0.543 
#>  9 viper     RFXAP  GSM2753335  0.488 
#> 10 viper     RFXAP  GSM2753336  1.32  
#> 11 viper     RFXAP  GSM2753337  1.93  
#> 12 viper     RFXAP  GSM2753338  1.93  
#> 13 viper     SMAD3  GSM2753335  0.176 
#> 14 viper     SMAD3  GSM2753336  0.0426
#> 15 viper     SMAD3  GSM2753337  0.219 
#> 16 viper     SMAD3  GSM2753338  0.142 
#> 17 viper     TFAP2A GSM2753335  0.722 
#> 18 viper     TFAP2A GSM2753336  0.582 
#> 19 viper     TFAP2A GSM2753337  0.462 
#> 20 viper     TFAP2A GSM2753338  0.330
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
