# decoupleR 1.2.0

## Changes
* Some method's names have been changed to make them easier to identify:
  * `pscira` now is called Weighted Sum (`wsum`).
  * `mean` now is called Weighted Mean (`wmean`).
  * `scira` now is called Univariate Linear Model (`ulm`).
  
* The column name for `tf` in the output tibbles has been changed to `source`.

* Updated documentation for all methods.

* Updated vignette and README.

* `decouple` function now accepts order mismatch between the list of methods and 
the list of methods's arguments.

* Moved benchmark branch to a separate repository as its own package: 
https://github.com/saezlab/decoupleRBench

## New features

* New methods added:
  * Fast Gene Set Enrichment Analysis (`fgsea`).
  * `AUCell`.
  * Univariate Decision Tree (`udt`).
  * Multivariate Decision Tree (`mdt`).
  * Multivariate Linear Model (`mlm`).

* New `decoupleR` manuscript repository: https://github.com/saezlab/decoupleR_manuscript

* New `consensus` score based on `RobustRankAggreg::aggregateRanks()` added when
running `decouple` with multiple methods. 

* Methods based on permutations or statistical tests now return also a p-value 
for the obtained score (`fgsea`, `mlm`, `ora`, `ulm`, `viper`, `wmean` and 
`wsum`).

* New error added when network edges are duplicated.

* New error added when the input matrix contains NAs or Infs. 

# decoupleR 1.1.0

## New features

All new features allow for **tidy selection**. Making it easier to evaluate
different types of data for the same method. For instance, you can specify the
columns to use as strings, integer position, symbol or expression.

### Methods

* New `decouple()` integrates the various member functions of the
  `decoupleR statistics` for centralized evaluation.
  
* New family `decoupleR statists` for shared documentation is made up of:
  * New `run_gsva()` incorporate a convinient wrapper for [GSVA::gsva()](https://rdrr.io/bioc/GSVA/man/gsva.html).
  * New `run_mean()` calculates both the unnormalized regulatory activity
    and the normalized (i.e. z-score) one based on an empirical distribution.
  * New `run_ora()` fisher exact test to calculate the regulatory activity.
  * New `run_pscira()` uses a logic equivalent to `run_mean()` with the
    difference that it does not accept a column of likelihood.
  * New `run_scira()` calculates the regulatory activity through the coefficient
    $\beta_1$ of an adjusted linear model.
  * New `run_viper()` incorporate a convinient wrapper for [viper::viper()](https://rdrr.io/bioc/viper/man/viper.html).

### Converters

* New functions family `convert_to_ variants` that allows the conversion
  of data to a standard format.
  * New `convert_to_()` return the entry without modification.
  * New `convert_to_gsva()` return a list of regulons suitable for [GSVA::gsva()](https://rdrr.io/bioc/GSVA/man/gsva.html).
  * New `convert_to_mean()` return a tibble with four columns:
    `tf`, `target`, `mor` and `likelihood`.
  * New `convert_to_ora()` returns a named list of regulons; tf with
    associated targets.
  * New `convert_to_pscira()` returns a tibble with three columns:
    `tf`, `target` and `mor`.
  * New `convert_to_scira()` returns a tibble with three columns:
    `tf`, `target` and `mor`.
  * New `convert_to_viper()` return a list of regulons suitable for
    [viper::viper()](https://rdrr.io/bioc/viper/man/viper.html)
