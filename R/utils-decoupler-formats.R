# File description --------------------------------------------------------
# The purpose of this file is to generate a documentation skeleton from
# which the methods in the package can inherit the corresponding attributes.
# Thus avoiding duplication of documentation while maintaining uniformity
# through the definition of functions in the package.

# mat --------------------------------------------------------------------

#' DecoupleR mat format
#'
#' @param mat Matrix to evaluate (e.g. expression matrix).
#'  Target nodes in rows and conditions in columns.
#'  `rownames(mat)` must have at least one intersection with the elements
#'  in `network` `.target` column.
#'
#' @keywords internal
#' @family decoupleR formats
#' @name .decoupler_mat_format
#' @aliases mat_format
NULL

# network -----------------------------------------------------------------

#' DecoupleR network format
#'
#' @description
#' A network passed to any `run_` method in the package must contain at
#' least two attributes: `.source` and `.target`. In addition,
#' the methods must map their corresponding metadata associated with their edges.
#'
#' @param network Tibble or dataframe with edges and it's associated metadata.
#' @param .source Column with source nodes.
#' @param .target Column with target nodes.
#' @param .mor Column with edge mode of regulation (i.e. mor).
#' @param .likelihood Column with edge likelihood.
#'
#' @details
#' * All the attributes to be mapped are prefixed by `.`
#' * The idea of using this type of mapping is to provide flexibility to
#'   different types of networks, be they regulatory, metabolic, or of any
#'   other type. This way, you should only consider having your network or
#'   networks in a long format and these can easily be manipulated by functions
#'   within the [tidyverse ecosystem](https://www.tidyverse.org/).
#'
#' @keywords internal
#' @family decoupleR formats
#' @name .decoupler_network_format
#' @aliases network_format
NULL


# benchmark input --------------------------------------------------------------
#' Benchmark input tibble containing the (experimental) design for each of the
#' benchmark runs corresponding to rows
#' @name input_tibble
#'
#' @details A tibble with locations, options, and filter options for
#' the desired benchmark setting
#'
#' @param set_name user-defined name of the set resource
#' @param bench_name user-defined name of the benchmark data
#' @param stats_list List of statistics to run
#' @param opts_list Named list containing the options for each stat. method
#' @param bexpr_loc benchmark expression data location (.rds format tibble)
#' @param bmeta_loc benchmark metadata location (.rds format tibble)
#' @param source_loc set source (e.g. network resource, gene ontology sets,
#'  kinase sets, etc.) location (.rds format tibble)
#' @param source_col name of the column with the source for the set source
#' @param target_col name of the column with the targets for the set source
#' @param filter_col name of the column by which we wish to filter
#' @param filter_crit criteria by which we wish to filter the `filter_col`
NULL
