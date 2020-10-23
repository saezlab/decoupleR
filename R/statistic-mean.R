#' SAFE mean
#'
#' Calculate the activity of a regulon through the conditions in the \code{mat}
#' matrix by calculating the mean over the expression of all genes.
#'
#' @param mat Evaluation matrix (e.g. expression matrix).
#'  Target nodes in rows and conditions in columns.
#' @param network Tibble or dataframe with edges and metadata.
#' @param .source Column with source nodes.
#' @param .target Column with target nodes.
#' @param .mor Column with edge mode of regulation (mor).
#' @param .likelihood Column with edge likelihood.
#' @param minsize How many output edges a source node must have to be included
#'  in the analysis?
#' @param times How many permutations to do?
#' @param seed A single value, interpreted as an integer, or NULL for random number generation.
#' @param sparse Should the matrices used for the calculation be sparse?
#' @param randomize_type How to randomize the expression matrix.
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the conditions. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item{\code{tf}}: {Source nodes of \code{network}.}
#'    \item{\code{condition}}: {Condition representing each column of \code{mat}.}
#'    \item{\code{score}}: {Regulatory activity (enrichment score).}
#'    \item{\code{statistic}}: {Indicates which method is associated with which score.}
#'    \item{\code{p_value}}: {p-value for the score of mean method.}
#'  }
#'
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
run_mean <- function(mat,
                     network,
                     .source = .data$tf,
                     .target = .data$target,
                     .mor = .data$mor,
                     .likelihood = .data$likelihood,
                     minsize = 1,
                     times = 2,
                     seed = 42,
                     sparse = TRUE,
                     randomize_type = "rows") {

  # Before to start ---------------------------------------------------------
  if (times < 2) {
    stop(str_interp("Parameter 'times' must be greater than or equal to 2, but ${times} was passed."))
  }

  network <- network %>%
    convert_to_mean({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Preprocessing -----------------------------------------------------------

  # Calculate the weights that will be used for the evaluation of the model
  network <- network %>%
    filter(.data$target %in% rownames(mat)) %>%
    add_count(.data$tf, name = "out_degree") %>%
    filter(.data$out_degree >= minsize) %>%
    add_count(.data$tf, wt = .data$likelihood, name = "contribution") %>%
    mutate(weight = .data$mor * .data$likelihood / .data$contribution) %>%
    select(-.data$out_degree, -.data$contribution)

  # Extract labels that will map to the expression and profile matrices
  tfs <- network %>%
    pull(.data$tf) %>%
    unique()

  shared_targets <- network %>%
    pull(.data$target) %>%
    unique()

  targets <- rownames(mat)
  conditions <- colnames(mat)

  # Extract matrix of weights
  weight_mat <- network %>%
    pivot_wider_profile(
      .data$tf,
      .data$target,
      .data$weight,
      to_matrix = TRUE,
      to_sparse = sparse,
      values_fill = 0
    )

  # Analysis ----------------------------------------------------------------
  .mean_analysis(mat, weight_mat, shared_targets, times, seed, randomize_type)
}

# Helper functions --------------------------------------------------------

#' Wrapper to execute run_mean() logic once finished preprocessing of data
#'
#' @inherit run_mean description
#'
#' @inheritParams run_mean
#' @param weight_mat Matrix that corresponds to the multiplication of the mor
#'  column with likelihood divided over the contribution.
#' @param shared_targets Target nodes that are shared between the
#'  \code{mat} and \code{network}.
#'
#' @inherit run_mean return
#'
#' @keywords internal
#' @noRd
.mean_analysis <- function(mat, weight_mat, shared_targets, times, seed, randomize_type) {
  # Thus, it is only necessary to define if we want
  # to evaluate a random model or not.
  mean_run <- partial(
    .mean_run,
    mat = mat,
    weight_mat = weight_mat,
    shared_targets = shared_targets,
    randomize_type = randomize_type
  )

  # Set a seed to ensure reproducible results
  set.seed(seed)
  # Run model for random data
  map_dfr(1:times, ~ mean_run(random = TRUE)) %>%
    group_by(.data$tf, .data$condition) %>%
    summarise(
      null_distribution = list(.data$value),
      null_mean = mean(.data$value),
      null_sd = stats::sd(.data$value),
      .groups = "drop"
    ) %>%
    # Run the true model and joined to random.
    left_join(y = mean_run(random = FALSE), by = c("tf", "condition")) %>%
    # Calculate scores
    mutate(
      z_score = (.data$value - .data$null_mean) / .data$null_sd,
      z_score = replace_na(.data$z_score, 0),
      p_value = map2_dbl(
        .x = .data$null_distribution,
        .y = .data$value,
        .f = ~ sum(abs(.x) > abs(.y)) / length(.x)
      )
    ) %>%
    # Reformat results
    select(-contains("null")) %>%
    rename(mean = .data$value, normalized_mean = .data$z_score) %>%
    pivot_longer(
      cols = c(.data$mean, .data$normalized_mean),
      names_to = "statistic",
      values_to = "score"
    ) %>%
    arrange(.data$statistic, .data$tf, .data$condition) %>%
    select(.data$tf, .data$condition, .data$score, .data$statistic, .data$p_value)
}

#' Wrapper to run mean one time
#'
#' @inheritParams .mean_analysis
#' @inherit .mean_evaluate_model return
#' @keywords internal
#' @noRd
.mean_run <- function(mat, weight_mat, shared_targets, random, randomize_type) {
  .mean_map_model_data(mat, shared_targets, random, randomize_type) %>%
    .mean_evaluate_model(weight_mat)
}

#' Collect a subset of data: random or not.
#'
#' If random is true, then it permutes the rows of the matrix
#' (i.e preserves the column relationships), otherwise it maintains
#' the original order of the data. Then it takes only those rows with
#' the values provided in \code{shared_targets}.
#'
#' @return Matrix with rows that match \code{shared_targets}.
#'
#' @inheritParams .mean_analysis
#' @keywords internal
#' @noRd
.mean_map_model_data <- function(mat, shared_targets, random, randomize_type) {
  if (random) {
    randomize_matrix(mat, randomize_type = randomize_type)[shared_targets, ]
  } else {
    mat[shared_targets, ]
  }
}

#' Evaluate model
#'
#' The evaluation model consists of evaluating the multiplication of the
#' weights by the factor of interest and comparing it against results
#' from permutations of the matrix of values of interest.
#'
#' @inheritParams .mean_analysis
#'
#' @return A dataframe with three columns:
#'  tf (source nodes), condition (colnames of mat) and value (score).
#'
#' @keywords internal
#' @noRd
.mean_evaluate_model <- function(mat, weight_mat) {
  (weight_mat %*% mat) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-.data$tf, names_to = "condition")
}
