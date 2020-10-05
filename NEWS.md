# decoupleR 0.0.1.9000

## New features

All new features allow for **tidy evaluation**. Making it easier to evaluate
different types of data for the same method.

### Methods

*  New `run_scira()` calculates the regulatory activity of a tf given the MoR of its targets.  
   **Regulatory activity** is equal to the t-value of the coefficient $\beta_{1}$
   of a linear regression of the expression profiles given a MoR.
* New `run_pscira()` calculates the regulatory activity of a tf given the MoR of its targets.  
  **Regulatory activity** is equal to the z-score resulting from comparing the
  expression value of the tf targets multiplied by their MoR with a 
  distribution of values associated with random targets.

### Converters

* New functions family `convert_to[method]()` that allows the conversion of data to a standard format.
  * `convert_to_scira()`
  * `convert_to_pscira()`

# decoupleR 0.0.0.9000
* Project kickoff
