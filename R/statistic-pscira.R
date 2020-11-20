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
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
run_pscira <- function(mat,
                       network,
                       .source = .data$tf,
                       .target = .data$target,
                       .mor = .data$mor,
                       sparse = TRUE,
                       times = 10,
                       seed = 42) {

  # Before to start ---------------------------------------------------------
  .start_time <- Sys.time()

  if (times < 2) {
    stop(str_interp("Parameter 'times' must be greater than or equal to 2, but ${times} was passed."))
  }

  # Preprocessing -----------------------------------------------------------

  # Convert to standard tibble: tf-target-mor.
  network <- network %>%
    convert_to_pscira({{ .source }}, {{ .target }}, {{ .mor }})

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
    pivot_wider_profile(.data$tf, .data$target, .data$mor, to_matrix = TRUE, to_sparse = sparse)

  # Convert to matrix to ensure that matrix multiplication works
  # in case mat is a labelled dataframe.
  mat <- as.matrix(mat)

  # Evaluate model ----------------------------------------------------------
  .pscira_analysis(mat, mor_mat, times, seed) %>%
    mutate(statistic_time = difftime(Sys.time(), .start_time))
}

# Helper functions --------------------------------------------------------

#' Wrapper to execute run_pscira() logic one finished preprocessing of data
#'
#' @inheritParams run_pscira
#' @param mor_mat Matrix that corresponds to the mor of the
#' target genes (columns) of a tf (rows).
#'
#' @inherit run_pscira return
#' @keywords intern
#' @noRd
.pscira_analysis <- function(mat, mor_mat, times, seed) {
  pscira_run <- partial(
    .f = .pscira_run,
    mat = mat,
    mor_mat = mor_mat
  )

  set.seed(seed)
  map_dfr(1:times, ~ pscira_run(random = TRUE)) %>%
    group_by(.data$tf, .data$condition) %>%
    summarise(.mean = mean(.data$value), .sd = sd(.data$value), .groups = "drop") %>%
    left_join(pscira_run(random = FALSE), by = c("tf", "condition")) %>%
    mutate(
      score = (.data$value - .data$.mean) / .data$.sd,
      score = replace_na(.data$score, 0)
    ) %>%
    transmute(statistic = "pscira", .data$tf, .data$condition, .data$score)
}

#'  Wrapper to perform mat %*% mor_mat
#'
#' @inheritParams .pscira_analysis
#' @param random Logical value that indicates whether the rows of the matrix
#' should be shuffled or not.
#'
#' @inherit .pscira_evaluate_model return
#' @keywords internal
#' @noRd
.pscira_run <- function(mat, mor_mat, random) {
  .pscira_map_model_data(mat, random) %>%
    .pscira_evaluate_model(mor_mat)
}

#' Map model data
#'
#' Build a data set with the necessary values to evaluate the model.
#'
#' @inheritParams .pscira_run
#'
#' @return origin nal/shuffled matrix
#' @keywords internal
#' @noRd
.pscira_map_model_data <- function(mat, random = FALSE) {
  if (random) {
    return(mat[sample(nrow(mat)), ])
  } else {
    mat
  }
}

#' Evaluate model
#'
#' Calculates the regulatory activity of all tfs with respect to its
#' associated profile for each condition.
#'
#' @inheritParams .pscira_run
#'
#' @return Tibble with tf regulatory activity for each tf-sample pair.
#' @keywords internal
#' @noRd
.pscira_evaluate_model <- function(mat, mor_mat) {
  (mor_mat %*% mat) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-.data$tf, names_to = "condition")
}
