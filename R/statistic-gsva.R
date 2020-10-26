#' GSVA wrapper
#'
#' This function is a convenient wrapper for the
#' \code{\link[=gsva]{GSVA::gsva()}} function.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param options A list of named options to pass to
#' \code{\link[=gsva]{GSVA::gsva()}}..
#' These options should not \code{include}, \code{expr} or \code{gset.idx.list}.
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the conditions. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item{\code{statistic}}: {Indicates which method is associated with which score.}
#'    \item{\code{tf}}: {Source nodes of \code{network}.}
#'    \item{\code{condition}}: {Condition representing each column of \code{mat}.}
#'    \item{\code{score}}: {Regulatory activity (enrichment score).}
#'  }
#' @export
run_gsva <- function(mat,
                     network,
                     .source = .data$tf,
                     .target = .data$target,
                     options = list()) {
  # Before to start ---------------------------------------------------------
  regulons <- network %>%
    convert_to_gsva({{ .source }}, {{ .target }})

  # Analysis ----------------------------------------------------------------
  do.call(
    what = GSVA::gsva,
    args = c(
      list(
        expr = mat,
        gset.idx.list = regulons
      ),
      options
    )
  ) %>%
    as.data.frame() %>%
    rownames_to_column(var = "tf") %>%
    pivot_longer(cols = -.data$tf, names_to = "condition", values_to = "score") %>%
    transmute(statistic = "gsva", .data$tf, .data$condition, .data$score)
}
