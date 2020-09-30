#' PSCIRA (Permutation Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculates TF activity according to \href{https://www.biorxiv.org/content/10.1101/553040v1.full.pdf}{SCIRA}.
#'
#' @inherit run_scira details
#'
#' @inheritParams run_scira
#' @param times Number of replications.
#' @param seed A single value, interpreted as an integer, or NULL.
#'
#' @inherit run_scira return
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr pull rename group_by summarise left_join mutate transmute
#' @importFrom purrr map_dfr
run_pscira <- function(emat,
                       regulons,
                       .source = .data$tf,
                       .target = .data$target,
                       .profile = .data$mor,
                       .sparse = TRUE,
                       times = 10,
                       seed = 42) {

  # Internal functions ------------------------------------------------------

  # Map model data
  #
  # Build a data set with the necessary values to evaluate the model.
  #
  # @param .emat Expression matrix.
  # @param random Logical value that indicates whether the rows of the matrix should be shuffled or not.
  #
  # @return Expression matrix.
  map_model_data <- function(.emat, random = FALSE) {
    if (random) {
      return(.emat[sample(nrow(.emat)), ])
    } else {
      .emat
    }
  }

  # Evaluate model
  #
  # Calculates the regulatory activity of all tfs with respect to its
  # associated profile for each condition.
  #
  # @param .emat Expression matrix.
  #
  # @return Tibble with tf regulatory activity for each tf-sample pair.
  evaluate_model <- function(.emat) {
    (profile_mat %*% .emat) %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("source") %>%
      pivot_longer(-.data$source, names_to = "condition")
  }

  # Preprocessing -----------------------------------------------------------

  # Extract labels that will map to the expression and profile matrices
  sources <- regulons %>%
    pull({{ .source }}) %>%
    unique()

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  profile_mat <- regulons %>%
    rename(source = {{ .source }}, target = {{ .target }}, profile = {{ .profile }}) %>%
    get_profile_of(
      sources = list(source = sources, target = rownames(emat)),
      values_fill = list(profile = 0)
    ) %>%
    pivot_wider_profile(.data$source, .data$target, .data$profile, to_matrix = TRUE, to_sparse = .sparse)

  # Convert to matrix to ensure that matrix multiplication works
  # in case emat is a labelled dataframe.
  emat <- as.matrix(emat)

  # Evaluate model ----------------------------------------------------------

  set.seed(seed)
  map_dfr(1:times, ~ map_model_data(emat, random = TRUE) %>%
            evaluate_model()) %>%
    group_by(.data$source, .data$condition) %>%
    summarise(.mean = mean(.data$value), .sd = stats::sd(.data$value), .groups = "drop") %>%
    left_join(evaluate_model(emat), by = c("source", "condition")) %>%
    mutate(score = (.data$value - .data$.mean) / .data$.sd) %>%
    transmute(.data$source, .data$condition, .data$score)
}
