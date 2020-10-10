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
#' @import dplyr
run_mean <- function(mat,
                     network,
                     .source = .data$tf,
                     .target = .data$target,
                     .mor = .data$mor,
                     .likelihood = .data$likelihood,
                     minsize = 1,
                     times = 2,
                     seed = 42,
                     sparse = TRUE) {

  # Before to start ---------------------------------------------------------
  if (times < 2) {
    stop(str_interp("Parameter 'times' must be greater than or equal to 2, but ${times} was passed."))
  }

  network <- network %>%
    convert_to_mean({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Preprocessing -----------------------------------------------------------

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
  samples <- colnames(mat)

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

  mean_run <- partial(
    .mean_run,
    .mat = mat,
    .weight_mat = weight_mat,
    .shared_targets = shared_targets
  )


  set.seed(seed)
  map_dfr(1:times, ~ mean_run(random = TRUE)) %>%
    group_by(.data$tf, .data$sample) %>%
    summarise(
      null_distribution = list(.data$value),
      null_mean = mean(.data$value),
      null_sd = stats::sd(.data$value),
      .groups = "drop"
    ) %>%
    left_join(y = mean_run(), by = c("tf", "sample")) %>%
    mutate(
      z_score = (.data$value - .data$null_mean) / .data$null_sd,
      z_score = replace_na(.data$z_score, 0),
      p_value = map2_dbl(
        .x = .data$null_distribution,
        .y = .data$value,
        .f = ~ sum(abs(.x) > abs(.y)) / length(.x)
      )
    ) %>%
    select(-contains("null")) %>%
    rename(mean = .data$value, normalized_mean = .data$z_score) %>%
    pivot_longer(
      cols = c(.data$mean, .data$normalized_mean),
      names_to = "statistic",
      values_to = "score"
    ) %>%
    arrange(.data$statistic, .data$tf, .data$sample) %>%
    select(.data$tf, .data$sample, .data$score, .data$statistic, .data$p_value)
}

# Helper functions --------------------------------------------------------

.mean_map_model_data <- function(.mat, .shared_targets, random = FALSE) {
  if (random) {
    .mat[sample(nrow(.mat)), ] %>%
      `rownames<-`(rownames(.mat)) %>%
      .[.shared_targets, ]
  } else {
    .mat[.shared_targets, ]
  }
}

.mean_evaluate_model <- function(.mat, .weight_mat) {
  (.weight_mat %*% .mat) %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-.data$tf, names_to = "sample")
}

.mean_run <- function(.mat, .weight_mat, .shared_targets, random = FALSE) {
  .mean_map_model_data(.mat, .shared_targets, random = random) %>%
    .mean_evaluate_model(.weight_mat)
}
