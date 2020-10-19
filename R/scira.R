#' SCIRA (Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculates TF activity according to \href{https://www.biorxiv.org/content/10.1101/553040v1.full.pdf}{SCIRA}.
#'
#' @details
#' Estimation of regulatory activity: A linear regression of the expression profile
#' is performed against the "target profile" of the given TF, where in the target
#' profile, any member of the regulon will be assigned a +1 for activated actions,
#' a -1 for inhibitory activations and, finally, a 0 for all genes that are not
#' a member of the TF regulon. TF activity is then defined as the t-statistic of
#' this linear regression.
#'
#' @param emat A named expression matrix. E.g Genes in rows and samples in columns.
#' @param regulons A data frame of regulons in table format.
#' @param .source Column name in regulons with the factors to calculate the enrichment scores.
#' @param .target Column name in regulons that relates the labels of the rows in the expression matrix.
#' @param .target_profile Column name in regulons that contains the target profile. E.g. MoR.
#' @param .sparse Logical value indicating if the generated profile matrix should be sparse.
#'
#' @return A long format tibble of the enrichment results for each set of genes
#'  across the samples. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item \code{source} Resources from the \code{.source} column of the \code{regulons} data frame.
#'    \item \code{condition} Samples representing each column of \code{emat}.
#'    \item \code{score} Regulatory activity of each resource for each condition.
#'  }
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats coef lm summary.lm
run_scira <- function(emat,
                      regulons,
                      .source = .data$tf,
                      .target = .data$target,
                      .target_profile = .data$mor,
                      .sparse = FALSE) {

  # Preprocessing -----------------------------------------------------------

  # Convert to standard tibble: tf-target-mor.
  regulons <- regulons %>%
    convert_to_scira({{ .source }}, {{ .target }}, {{ .target_profile }}, clean = TRUE)

  # Extract labels that will map to the expression and profile matrices
  tfs <- regulons %>%
    pull(.data$tf) %>%
    unique()

  conditions <- colnames(emat)

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  profile_mat <- regulons %>%
    get_profile_of(
      sources = list(tf = tfs, target = rownames(emat)),
      values_fill = list(mor = 0)
    ) %>%
    pivot_wider_profile(.data$tf, .data$target, .data$mor, to_matrix = TRUE, to_sparse = .sparse)

  # Model evaluation --------------------------------------------------------
  .scira_analysis(emat, profile_mat)
}

# Helper functions ------------------------------------------------------

#' Wrapper to execute run_scira() logic one finished preprocessing of data
#'
#' @inheritParams run_mean
#' @param target_profile_mat
#'
#' @inherit run_mean return
#' @keywords intern
#' @noRd
.scira_analysis <- function(emat, target_profile_mat) {

  # Allocate the space for all combinations of sources and conditions
  # and evaluate the proposed model.
  lift_dl(expand_grid)(list(tf = rownames(target_profile_mat), condition = colnames(emat))) %>%
    rowwise() %>%
    mutate(score = .scira_run(.data$tf, .data$condition, emat, target_profile_mat))
}

#' Wrapper to run scira one tf per sample at time
#'
#' @inheritParams .scira_analysis
#' @param source Current `source` to evaluate.
#' @param condition Current `condition` to evaluate
#'
#' @inherit .scira_evaluate_model return
#'
#' @keywords internal
#' @noRd
.scira_run <- function(source, condition, emat, target_profile_mat) {
  .scira_map_model_data(source, condition, emat, target_profile_mat) %>%
    .scira_evaluate_model()
}

#' Map model data
#'
#' Build a data set with the necessary values to evaluate the model.
#'
#' @inheritParams .scira_run
#'
#' @return Tibble with the data to use to evaluate the model.
#' @keywords internal
#' @noRd
.scira_map_model_data <- function(source, condition, emat, target_profile_mat) {
  tibble(expression = emat[, condition], mor = target_profile_mat[source, ])
}

#' Evaluate model
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @param data A data frame, data frame extension (e.g. a tibble) with two
#'  columns; `expression` and `mor` to perform lineal regression.
#'
#' @return t-value corresponding to beta 1 parameter of linear regression.
#' @keywords internal
#' @noRd
.scira_evaluate_model <- function(data) {
  coefficients_values <- lm(expression ~ mor, data = data) %>%
    summary.lm() %>%
    coef()

  t_values <- coefficients_values[, "t value"]
  pluck(t_values, 2, .default = NA)
}
