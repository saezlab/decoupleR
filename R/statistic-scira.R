#' SCIRA (Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculates TF activity according to [SCIRA](https://www.biorxiv.org/content/10.1101/553040v1.full.pdf).
#'
#' @details
#' Estimation of regulatory activity: A linear regression of the expression profile
#' is performed against the "target profile" of the given TF, where in the target
#' profile, any member of the regulon will be assigned a +1 for activated actions,
#' a -1 for inhibitory activations and, finally, a 0 for all genes that are not
#' a member of the TF regulon. TF activity is then defined as the t-statistic of
#' this linear regression.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Logical value indicating if the generated profile matrix should be sparse.
#' @param fast Logical value indicating if the lineal model must be calculated
#' with [speedglm::speedlm.fit()] or with base
#' [stats::lm()].
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the conditions. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#'  5. `statistic_time`: Internal execution time indicator.
#' @family decoupleR statistic
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats coef lm summary.lm
#' @importFrom speedglm speedlm.fit
run_scira <- function(mat,
                      network,
                      .source = .data$tf,
                      .target = .data$target,
                      .mor = .data$mor,
                      sparse = FALSE,
                      fast = TRUE) {

  # Preprocessing -----------------------------------------------------------
  .start_time <- Sys.time()

  # Convert to standard tibble: tf-target-mor.
  network <- network %>%
    convert_to_scira({{ .source }}, {{ .target }}, {{ .mor }})

  # Extract labels that will map to the expression and profile matrices
  tfs <- network %>%
    pull(.data$tf) %>%
    unique()

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  mor_mat <- network %>%
    get_profile_of(
      sources = list(tf = tfs, target = rownames(mat)),
      values_fill = list(mor = 0)
    ) %>%
    pivot_wider_profile(.data$target, .data$tf, .data$mor, to_matrix = TRUE, to_sparse = sparse)

  # Model evaluation --------------------------------------------------------
  .scira_analysis(mat, mor_mat, fast) %>%
    mutate(statistic_time = difftime(Sys.time(), .start_time))
}

# Helper functions ------------------------------------------------------

#' Wrapper to execute run_scira() logic one finished preprocessing of data
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @inheritParams run_scira
#' @param mor_mat
#'
#' @inherit run_scira return
#' @keywords intern
#' @noRd
.scira_analysis <- function(mat, mor_mat, fast) {
  scira_evaluate_model <- partial(
    .f = .scira_evaluate_model,
    mat = mat,
    mor_mat = mor_mat,
    fast = fast
  )

  # Allocate the space for all combinations of sources and conditions
  # and evaluate the proposed model.
  lift_dl(expand_grid)(list(tf = colnames(mor_mat), condition = colnames(mat))) %>%
    rowwise(.data$tf, .data$condition) %>%
    summarise(
      score = scira_evaluate_model(.data$tf, .data$condition),
      .groups = "drop"
    ) %>%
    transmute(statistic = "scira", .data$tf, .data$condition, .data$score)
}

#' Wrapper to run scira one tf (source) per sample (condition) at time
#'
#' @keywords internal
#' @noRd
.scira_evaluate_model <- function(source, condition, mat, mor_mat, fast) {
  if (fast) {
    speedlm.fit(
      y = mat[, condition],
      X = cbind(1, mor_mat[, source])
    ) %>%
      summary() %>%
      pluck("coefficients", "t", 2, .default = NA)
  } else {
    coefficients_values <- lm(mat[, condition] ~ mor_mat[, source]) %>%
      summary() %>%
      coef()

    t_values <- coefficients_values[, "t value"]
    pluck(t_values, 2, .default = NA)
  }
}
