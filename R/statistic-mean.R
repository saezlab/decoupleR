#' Weighted mean
#'
#' Calculate the activity of all regulons in `network` through the conditions in
#' the `mat` matrix by calculating the mean over the expression of all genes.
#'
#' @details
#'  `run_mean()` calculates the activity score, but in addition, it takes
#'  advantage of the permutations used to calculate the `p-value`, to provide
#'  the normalized activity score. This is represented in the `statistic` column
#'  which will contain two values for each call to `run_mean()`; __mean__ and
#'  __normalized_mean__.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param times How many permutations to do?
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
#' @param sparse Should the matrices used for the calculation be sparse?
#' @param randomize_type How to randomize the expression matrix.
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#'  5. `p_value`: p-value for the score of mean method.
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats sd
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_mean(mat, network, tf, target, mor, likelihood)
run_mean <- function(mat,
                     network,
                     .source = .data$tf,
                     .target = .data$target,
                     .mor = .data$mor,
                     .likelihood = .data$likelihood,
                     times = 2,
                     seed = 42,
                     sparse = TRUE,
                     randomize_type = "rows") {
    # Before to start ---------------------------------------------------------
    if (times < 2) {
        stop(stringr::str_interp("Parameter 'times' must be greater than or equal to 2, but ${times} was passed."))
    }

    network <- network %>%
        convert_to_mean({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

    # Preprocessing -----------------------------------------------------------

    # Calculate the weights that will be used for the evaluation of the model
    network <- network %>%
        filter(.data$target %in% rownames(mat)) %>%
        .mean_calculate_weight()

    # Extract labels that will map to the expression and profile matrices
    tfs <- unique(network[["tf"]])

    shared_targets <- unique(network[["target"]])

    targets <- rownames(mat)
    conditions <- colnames(mat)

    # Extract matrix of weights
    weight_mat <- network %>%
        pivot_wider_profile(
            id_cols = .data$tf,
            names_from = .data$target,
            values_from = .data$weight,
            to_matrix = TRUE,
            to_sparse = sparse,
            values_fill = 0
        )

    # Analysis ----------------------------------------------------------------
    withr::with_seed(seed, {
        .mean_analysis(mat, weight_mat, shared_targets, times, seed, randomize_type)
    })

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
#'  `mat` and `network`.
#'
#' @inherit run_mean return
#'
#' @keywords internal
#' @noRd
.mean_analysis <- function(mat, weight_mat, shared_targets, times, randomize_type) {
    # Thus, it is only necessary to define if we want
    # to evaluate a random model or not.
    mean_run <- partial(
        .mean_run,
        mat = mat,
        weight_mat = weight_mat,
        shared_targets = shared_targets,
        randomize_type = randomize_type
    )

    # Run model for random data
    map_dfr(seq_len(times), ~ mean_run(random = TRUE)) %>%
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
        select(.data$statistic, .data$tf, .data$condition, .data$score, .data$p_value)
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

#' Calculate mean weight
#'
#' @inheritParams .mean_analysis
#' @keywords internal
#' @noRd
.mean_calculate_weight <- function(network) {
    network %>%
        add_count(.data$tf, name = "contribution") %>%
        transmute(
            .data$tf,
            .data$target,
            weight = .data$mor * .data$likelihood / .data$contribution
        )
}

#' Collect a subset of data: random or not.
#'
#' If random is true, then it permutes the rows of the matrix
#' (i.e preserves the column relationships), otherwise it maintains
#' the original order of the data. Then it takes only those rows with
#' the values provided in `shared_targets`.
#'
#' @return Matrix with rows that match `shared_targets`.
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
