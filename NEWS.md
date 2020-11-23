# decoupleR 0.0.1.9000

## New features

All new features allow for **tidy selection**. Making it easier to evaluate
different types of data for the same method. For instance, you can specify the
columns to use as strings, integer position, symbol or expression.

### Methods

* The new `decouple()` integrates the various member functions of the
  `decoupleR statistics` for centralized evaluation.
  
* New family `decoupleR statists` for shared documentation is made up of:
  * New `run_gsva()`
  * New `run_mean()`
  * New `run_pscira()`
  * New `run_scira()`
  * New `run_viper()`

### Converters

* New functions family `convert_to_ variants` that allows the conversion of data to a standard format.
  * New `convert_to_()` return the entry without modification.
  * New `convert_to_gsva()` return a list of regulons suitable for GSVA::gsva().
  * New `convert_to_mean()` return a tibble with four columns: `tf`, `target`, `mor` and `likelihood`.
  * New `convert_to_pscira()` returns a tibble with three columns: `tf`, `target` and `mor`.
  * New `convert_to_scira()` returns a tibble with three columns: `tf`, `target` and `mor`.
  * New `convert_to_viper()` return a list of regulons suitable for viper::viper()

# decoupleR 0.0.0.9000
* Project kickoff
