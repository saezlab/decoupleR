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
#'  \code{rownames(mat)} must have at least one intersection with the elements
#'  in \code{network} \code{.target} column.
#'
#' @keywords internal
#' @name .decoupler_mat_format
NULL

# network -----------------------------------------------------------------

#' DecoupleR network format
#'
#' @description
#' A network passed to any \code{run_} method in the package must contain at
#' least two attributes: \code{.source} and \code{.target}. In addition,
#' the methods must map their corresponding metadata associated with their edges.
#'
#' @param network Tibble or dataframe with edges and it's associated metadata.
#' @param .source Column with source nodes.
#' @param .target Column with target nodes.
#' @param .mor Column with edge mode of regulation (i.e. mor).
#' @param .likelihood Column with edge likelihood.
#' @keywords internal
#' @name .decoupler_network_format
NULL
