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
#' @param mat A named expression matrix. E.g Genes in rows and samples in columns.
#' @param network Tibble or dataframe with edges and associated metadata.
#' @param .source Column with source nodes.
#' @param .target Column with target nodes.
#' @param .mor Column with edge mode of regulation (mor).
#' @param .sparse Logical value indicating if the generated profile matrix should be sparse.
#'
#' @return A long format tibble of the enrichment results for each set of genes
#'  across the samples. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item \code{source} Resources from the \code{.source} column of the \code{network} data frame.
#'    \item \code{condition} Samples representing each column of \code{mat}.
#'    \item \code{score} Regulatory activity of each resource for each condition.
#'  }
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats coef lm summary.lm
run_scira <- function(mat,
                      network,
                      .source = .data$tf,
                      .target = .data$target,
                      .mor = .data$mor,
                      .sparse = FALSE) {

  # Preprocessing -----------------------------------------------------------

  # Convert to standard tibble: tf-target-mor.
  network <- network %>%
    convert_to_scira({{ .source }}, {{ .target }}, {{ .mor }})

  # Extract labels that will map to the expression and profile matrices
  tfs <- network %>%
    pull(.data$tf) %>%
    unique()

  conditions <- colnames(mat)

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  profile_mat <- network %>%
    get_profile_of(
      sources = list(tf = tfs, target = rownames(mat)),
      values_fill = list(mor = 0)
    ) %>%
    pivot_wider_profile(.data$tf, .data$target, .data$mor, to_matrix = TRUE, to_sparse = .sparse)

  # Model evaluation --------------------------------------------------------
  .scira_analysis(mat, profile_mat)
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
.scira_analysis <- function(mat, mor_mat) {

  # Allocate the space for all combinations of sources and conditions
  # and evaluate the proposed model.
  lift_dl(expand_grid)(list(tf = rownames(mor_mat), condition = colnames(mat))) %>%
    rowwise() %>%
    mutate(
      score = .scira_run(.data$tf, .data$condition, mat, mor_mat),
      statistic = "scira"
    ) %>%
    ungroup()
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
.scira_run <- function(source, condition, mat, mor_mat) {
  .scira_map_model_data(source, condition, mat, mor_mat) %>%
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
.scira_map_model_data <- function(source, condition, mat, mor_mat) {
  tibble(expression = mat[, condition], mor = mor_mat[source, ])
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
