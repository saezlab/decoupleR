# Contributing to decoupleR

This outlines how to propose a change to decoupleR. 
For more detailed info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib). 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("saezlab/decoupleR", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  
   
## Add a new statistic

Read carefully the instructions dedicated to the pull request section above and consider the following specific steps:

*   Open an [issue](https://github.com/saezlab/decoupleR/issues) stating which statistic you would like to add.
    Assign it to yourself and label it with the badge `enhancement`; check before if not someone else is
    already working on implementing this statistic.  

*   Define the name of the statistic, for now, assume you want to define the statistic `foo`.
    
*   Create a file called `R/statistic-foo.R` and define
    `run_foo(mat, network, .source, .target, ...)` or
    `run_foo(mat, network, .source, .target, ..., options = list())`,
    depending if the statistic is designed from zero or is a wrapper, respectively;
    where `...` denote any extra arguments needed for the statistic.
    If the algorithm requires helper functions, add them here too.  
    **Return** of any function of the `run_` family must always be a tidy tibble with the following columns:  
    1. **statistic**: Indicates which methods is associated with which score.
    2. **tf**: Source nodes of `network`.
    3. **condition**: Conditions representing each column of `mat`.
    4. **score**: Regulatory activity (enrichment score).
    5. **statistic_time**: Internal execution time indicator.
    6. **...**: If the algorithm requires it, add the generated metadata columns. For instance, p-value.
      
    **Notes:**
    *   Check
        [utils-decoupler-formats.R](https://github.com/saezlab/decoupleR/blob/documentation/R/utils-decoupler-formats.R)
        to understand the arguments and specific conventions.
    *   For examples of implementations, please check
        [statistic-scira.R](https://github.com/saezlab/decoupleR/blob/documentation/R/statistic-scira.R) or
        [statistic-viper.R](https://github.com/saezlab/decoupleR/blob/documentation/R/statistic-viper.R)
        for design from scratch or wrapper, respectively.  
*   Inside
    [utils-dataset-converters.R](https://github.com/saezlab/decoupleR/blob/documentation/R/utils-dataset-converters.R)
    define `convert_to_foo(dataset, .source, .target, ...)`, where `...` indicate
    any extra *edge attribute* that needs to be mapped into the network.  
    **Goals:**  
    *   Transform the data frame that represents
        the network to a standardized format, so that, in the main function `run_foo()`,
        we can use and manipulate the columns to our liking without worrying about the mapping.  
    *   It allows a uniform mapping through the statistics contained in the package.
        What ensures rapid integration into the DecoupleR statistics ecosystem.  
      
    **Template** with the essential structure that the function must have.  
    ```r
    # foo ---------------------------------------------------------------------

    #' @rdname convert_to_
    #'
    #' @inheritParams run_foo
    #'
    #' @export
    #' @family convert_to_ variants
    convert_to_foo <- function(dataset, .source, .target) {
      .check_quos_status({{ .source }}, {{ .target }}, .dots_names = c(".source", ".target"))
    
      dataset %>%
        convert_f_defaults(
          tf = {{ .source }},
          target = {{ .target }}
        )
        # If extra conversion is needed add it here...
    }
    ```
*   Add the unit test for `run_foo()` only, once they pass incorporate `run_foo()`
    inside of the main `decouple()` and add it also to the `decouple()` unit testing.
*   If all goes well, you are ready. Start the pull request process.

## Code of Conduct

Please note that the decoupleR project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
