#' SAFE mean
#'
#' Calculate the activity of a regulon through the samples in the \code{mat}
#' matrix by calculating the mean over the expression of all genes.
#'
#' @param mat Evaluation matrix (e.g. expression matrix).
#'  Target nodes in rows and samples in columns.
#' @param network Tibble or dataframe with edges and metadata.
#' @param .source Column with source nodes.
#' @param .target Column with target nodes.
#' @param .mor Column with edge mode of regulation (mor).
#' @param .likelihood Column with edge likelihood.
#' @param minsize How many output edges a source node must have to be included
#'  in the analysis?
#' @param mixed If TRUE permutations are performed only within the positive and negative values.
#' @param times How many permutations to do?
#' @param seed A single value, interpreted as an integer, or NULL for random number generation.
#' @param sparse Should the matrices used for the calculation be sparse?
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item{\code{tf}}: {Source nodes of \code{network}}.
#'    \item{\code{sample}}: {Samples representing each column of \code{mat}}.
#'    \item{\code{score}}: {Regulatory activity (enrichment score)}.
#'  }
#'
#' @export
run_mean <- function(mat,
                     network,
                     .source = .data$tf,
                     .target = .data$target,
                     .mor = .data$mor,
                     .likelihood = .data$likelihood,
                     minsize = 1,
                     mixed = FALSE,
                     times = 2,
                     seed = 42,
                     sparse = FALSE) {

  # Before to start ---------------------------------------------------------
  if (times < 2) {
    stop(str_interp("Parameter 'times' must be greater than or equal to 2, but ${times} was passed."))
  }

  network <- network %>%
    convert_to_mean({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Preprocessing -----------------------------------------------------------


  # Analysis ----------------------------------------------------------------
}
