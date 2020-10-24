#' VIPER wrapper
#'
#' This function is a convenient wrapper for the
#' \code{\link[=viper]{viper::viper()}} function.
#'
#' @param emat An expression matrix with genes (HGNC symbol) in rows and samples
#'  in columns.
#' @param genesets A data frame of gene sets. The structure is dependent on the
#' gene set resource.
#' @param .source Column with source nodes.
#' @param .target Column with target nodes.
#' @param .mor Column with edge mode of regulation (mor).
#' @param .likelihood Column with edge likelihood.
#' @param options A list of named options to pass to
#' \code{\link[=viper]{viper::viper()}} such as \code{minsize} or \code{method}.
#'  These options should not \code{include}, \code{eset} or \code{regulon}.
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the conditions. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item{\code{statistic}}: {Indicates which method is associated with which score.}
#'    \item{\code{tf}}: {Source nodes of \code{network}.}
#'    \item{\code{condition}}: {Conditions representing each column of \code{mat}.}
#'    \item{\code{score}}: {Regulatory activity (enrichment score).}
#'  }
#'
#' @export
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import purrr
#' @import tidyr
#' @import viper
run_viper <- function(emat,
                      genesets,
                      .source = .data$tf,
                      .target = .data$target,
                      .mor = .data$mor,
                      .likelihood = .data$likelihood,
                      options = list()) {
  genesets <- genesets %>%
    convert_to_viper({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  do.call(
    viper,
    c(
      list(
        eset = emat,
        regulon = make_viper_genesets(genesets)
      ),
      options
    )
  ) %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-.data$tf, names_to = "condition", values_to = "score") %>%
    add_column(statistic = "viper", .before = 1)
}

#' Make gene sets for VIPER
#'
#' This function convert gene sets in a table format to the format required by
#' the \code{\link[=viper]{viper::viper()}} function.
#'
#' @param genesets A dataframe of gene sets that must contain the
#' columns \code{geneset}, \code{gene}, \code{mor} and \code{likelihood}.
#'
#' @return Gene sets in the \code{\link[=viper]{viper::viper()}} format.
#'
#' @keywords internal
#' @importFrom stats setNames
make_viper_genesets <- function(genesets) {
  genesets %>%
    split(.$geneset) %>%
    map(function(gs) {
      targets <- setNames(gs$mor, gs$gene)
      likelihood <- gs$likelihood
      list(tfmode = targets, likelihood = likelihood)
    })
}
