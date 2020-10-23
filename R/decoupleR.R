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
#' @param .options A list of argument-lists the same length as \code{statistics} (or length 1).
#'  The default argument, list(NULL), will be recycled to the same length as \code{statistics},
#'  and will call each function with no arguments (apart from \code{mat},
#'  \code{network}, \code{.source} and, \code{.target}).
#' @param statistics Statistical methods to be coupled.
#' @param run_ids What identifier to associate with each model evaluated?
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item{\code{run_id}}: {Indicates which statistic run is associeted to each observation.}
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
                     .options = list(NULL),
                     statistics,
                     run_ids = NULL) {

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

  # If requested add user identifiers to runs.
  if (!is_null(run_ids)) {
    statistics <- set_names(statistics, run_ids)
  }

  # Check options -----------------------------------------------------------
  if (is_empty(.options)) {
    .options <- list(NULL)
  }

  # Evaluate statistics -----------------------------------------------------

  # For the moment this will only ensure that the parameters passed
  # to decoupleR are the same when invoking the functions.
  invoke_map(
    .f = statistics,
    .x = .options,
    mat = mat,
    network = network,
    .source = enquo(.source),
    .target = enquo(.target)
  ) %>%
    bind_rows(.id = "run_id")
}

# Helpers -----------------------------------------------------------------
