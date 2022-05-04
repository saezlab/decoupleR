# decouplerR 2.2.

## Changes
* Changed example `mat` and `net` to toy examples.

* Changed test data to toy data.

## Bugfixes
* `ora` now selects correctly the top and bottom genes for p-value estimation.

# decoupleR 2.1.

## Changes
* `likelihood` param is deprecated, from now on, weights (positive or negative) 
  should go to the `mor` column of `network`. Methods will still run if 
  `likelihood` is specified, however they will be set to 1.

* Added `minsize` argument to all methods, set to 5 by default. Sources 
containing less than this value of targets in the input mat will be removed 
from the  calculations.

* Changed default behavior of the `decouple` function. Now if no methods are 
specified in the `statistics` argument, the function will only run the top 
performers in our benchmark (`mlm`, `ulm` and `wsum`). To run all methods like
before, set `statistics` to 'all'. Moreover, the argument `consensus_stats` has 
been added to filter statistics for the calculation of the `consensus` score. 
By default it only uses `mlm`, `ulm` and `norm_wsum`, or if `statistics`=='all'
all methods returned after running `decouple`.

* `viper` method:
    * Now properly handles weights in `mor` by normalizing them to -1 and +1.

* `ulm`/`mlm`/`udt`/`mdt` methods:
    * Changed how they processed the input network. Before the model 
    matrix only contained the intersection of features between mat and 
    network's targets, now it incorporates all features coming from mat 
    ensuring a more robust prediction. Prediction values may change slightly 
    from older versions. 
    * Deprecated `sparse` argument. 
    
* `ora` method:
    * Now takes top 5% features as default input instead of 300 up and bottom 
    features.
    * Added seed to randomly break ties
    
* `consensus` method: 
    * No longer based on `RobustRankAggreg`. Now the consensus score is the mean of the
    activities obtained after a double tailed z-score transformation.

* Discarded `filter_regulons` function.

* Moved major dependencies to Suggest to reduce the number of dependencies 
needed.

* Updated README by adding:
    * Kinase inference example
    * Graphical abstract
    * Manuscript and citation
    * New vignette style
    
* Updated documentation for all methods.

## New features
* Added wrappers to easily query `Omnipath`, one of the largest data-bases 
collecting prior-knowledge resources. Added these functions:
    * `show_resources`: shows available resources inside `Omnipath`.
    * `get_resource`: gets any resource from `Omnipath`.
    * `get_dorothea`: gets the DoRothEA gene regulatory network for 
    transcription factor (TF) activity estimation. Note: this version is 
    slightly different from the one in the package `dorothea` since it contains 
    new edges and TFs and also weights the interactions by confidence levels.
    * `get_progeny`: gets the PROGENy model for pathway activity estimation.

* Added `show_methods` function, it shows how many statistics are currently 
available.

* Added `check_corr` function, it shows how correlated regulators in a network 
are. It can be used to check for co-linearity for `mlm` and `mdt`. 

* Added new error for `mlm` when co-variables are co-linear (regulators are too 
correlated to fit a model).

## Bugfixes
* `wmean` and `wsum` now return the correct empirical p-values.

* `ulm`, `mlm`, `mdt` and `udt` now accept matrices with one column as input. 

* Results from `ulm` and `mlm` now correctly return un-grouped.

* Methods correctly run when `mat` has no column names.

# decoupleR 2.0

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

* New statistic `corr_wmean` inside `wmean`. 

* Methods based on permutations or statistical tests now return also a p-value 
for the obtained score (`fgsea`, `mlm`, `ora`, `ulm`, `viper`, `wmean` and 
`wsum`).

* New error added when network edges are duplicated.

* New error added when the input matrix contains NAs or Infs. 

# decoupleR 1.1

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
