# decoupleR 0.99.0

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
