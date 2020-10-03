#' PSCIRA (Permutation Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculate the regulatory activity of each tf by multiplying the expression
#' values of its objectives with their corresponding associated profiles for
#' each given condition.The result is equal to the z-score of the found value
#' compared to its null distribution.
#'
#' @inherit run_scira details
#'
#' @inheritParams run_scira
#' @param times Number of replications.
#' @param seed A single value, interpreted as an integer, or NULL.
#'
#' @inherit run_scira return
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
run_pscira <- function(emat,
                       regulons,
                       .source = .data$tf,
                       .target = .data$target,
                       .target_profile = .data$mor,
                       .sparse = TRUE,
                       times = 10,
                       seed = 42) {

  # Before to start ---------------------------------------------------------
  if (times < 2) {
    stop(str_interp("Parameter 'times' must be greater than or equal to 2, but ${times} was passed."))
  }

  # Preprocessing -----------------------------------------------------------

  # Convert to standard tibble: tf-target-mor.
  regulons <- regulons %>%
    convert_to_scira({{ .source }}, {{ .target }}, {{ .target_profile }}, clean = TRUE)

  # Extract labels that will map to the expression and profile matrices
  tfs <- regulons %>%
    pull(.data$tf) %>%
    unique()

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  target_profile_mat <- regulons %>%
    get_profile_of(
      sources = list(tf = tfs, target = rownames(emat)),
      values_fill = list(mor = 0)
    ) %>%
    pivot_wider_profile(.data$tf, .data$target, .data$mor, to_matrix = TRUE, to_sparse = .sparse)

  # Convert to matrix to ensure that matrix multiplication works
  # in case emat is a labelled dataframe.
  emat <- as.matrix(emat)

  # Evaluate model ----------------------------------------------------------

  set.seed(seed)
  map_dfr(1:times, ~ .pscira_map_model_data(emat, random = TRUE) %>%
    .pscira_evaluate_model(target_profile_mat)) %>%
    group_by(.data$tf, .data$condition) %>%
    summarise(.mean = mean(.data$value), .sd = sd(.data$value), .groups = "drop") %>%
    left_join(.pscira_evaluate_model(emat, target_profile_mat), by = c("tf", "condition")) %>%
    mutate(score = (.data$value - .data$.mean) / .data$.sd,
           score = replace_na(.data$score, 0)) %>%
    transmute(.data$tf, .data$condition, .data$score)
}

# Helper functions --------------------------------------------------------

#' Map model data
#'
#' Build a data set with the necessary values to evaluate the model.
#'
#' @param .emat Expression matrix.
#' @param random Logical value that indicates whether the rows of the matrix should be shuffled or not.
#'
#' @return Expression matrix.
#' @keywords internal
#' @noRd
.pscira_map_model_data <- function(.emat, random = FALSE) {
  if (random) {
    return(.emat[sample(nrow(.emat)), ])
  } else {
    .emat
  }
}

#' Evaluate model
#'
#' Calculates the regulatory activity of all tfs with respect to its
#' associated profile for each condition.
#'
#' @param .emat Expression matrix.
#' @param .target_profile_mat Matrix with specified profile.
#'
#' @return Tibble with tf regulatory activity for each tf-sample pair.
#' @keywords internal
#' @noRd
.pscira_evaluate_model <- function(.emat, .target_profile_mat) {
  (.target_profile_mat %*% .emat) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-.data$tf, names_to = "condition")
}
