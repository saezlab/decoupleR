#' decouple
#'
#' Calculate the TF activity per sample out of a gene expression matrix by
#' coupling a regulon network with a variety of statistics.
#'
#' @param mat Evaluation matrix (e.g. expression matrix).
#'  Target nodes in rows and samples in columns.
#' @param network Tibble or dataframe with edges and metadata.
#' @param .source Column in network with source nodes.
#' @param .target Column in network with target nodes.
#' @param .options Named list with edge attributes to use in the statistics.
#' @param statistics Statistical methods to be coupled.
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item{\code{tf}}: {Source nodes of \code{network}.}
#'    \item{\code{sample}}: {Samples representing each column of \code{mat}.}
#'    \item{\code{score}}: {Regulatory activity (enrichment score).}
#'    \item{\code{statistic}}: {Indicates which method is associated with which score.}
#'    \item{\code{metadata}}: {Metadata corresponding to the statistic collapsed to a string.}
#'  }
#' @export
#' @import purrr
decouple <- function(mat,
                     network,
                     .source,
                     .target,
                     .options = list(),
                     statistics) {

  # Match statistics to couple ----------------------------------------------

  # TODO this probably has to be changed to call the corresponding internal
  # functions of each statistic.
  available_statistics <- list(
    mean = run_mean,
    scira = run_scira,
    pscira = run_pscira,
    viper = run_viper
  )

  statistics <- statistics %>%
    match.arg(names(available_statistics), several.ok = TRUE) %>%
    available_statistics[.]

  # Check options -----------------------------------------------------------

  # For the moment this will only ensure that the parameters passed
  # to decoupleR are the same when invoking the functions.
  # TODO add function to check extra parameters.
  # TODO add function that allows list of list modification to ensure
  # shared parameters across statistics.
  .options <- list_modify(
    .x = .options,
    mat = mat,
    network = network,
    .source = ensym(.source),
    .target = ensym(.target)
  )

  # Evaluate statistics -----------------------------------------------------
  invoke_map(statistics, list(.options))
}

# Helpers -----------------------------------------------------------------
