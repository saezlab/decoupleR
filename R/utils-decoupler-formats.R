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
#' @param .likelihood Deprecated argument. Now it will always be set to 1.
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
