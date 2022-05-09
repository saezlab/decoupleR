#' Weighted Mean (WMEAN)
#'
#' @description
#' Calculates regulatory activities using WMEAN.
#'
#' @details
#' WMEAN infers regulator activities by first multiplying each target feature by
#' its associated weight which then are summed to an enrichment score
#' `wmean`. Furthermore, permutations of random target features can
#' be performed to obtain a null distribution that can be used to compute a
#' z-score `norm_wmean`, or a corrected estimate `corr_wmean` by multiplying
#' `wmean` by the minus log10 of the obtained empirical p-value.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param times How many permutations to do?
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
#' @param sparse Should the matrices used for the calculation be sparse?
#' @param randomize_type How to randomize the expression matrix.
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#'  5. `p_value`: p-value for the score of the method.
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
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_wmean(mat, net, minsize=0)
run_wmean <- function(mat,
                     network,
                     .source = .data$source,
                     .target = .data$target,
                     .mor = .data$mor,
                     .likelihood = .data$likelihood,
                     times = 100,
                     seed = 42,
                     sparse = TRUE,
                     randomize_type = "rows",
                     minsize = 5
                     ) {
    # Before to start ---------------------------------------------------------
    if (times < 2) {
        rlang::abort(message = stringr::str_glue("Parameter 'times' must be greater than or equal to 2, but {times} was passed."))
    }

    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    network <- network %>%
        rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
    network <- filt_minsize(rownames(mat), network, minsize)

    # Preprocessing -----------------------------------------------------------

    # Calculate the weights that will be used for the evaluation of the model
    network <- network %>%
        filter(.data$target %in% rownames(mat)) %>%
        .wmean_calculate_weight()

    # Extract labels that will map to the expression and profile matrices
    shared_targets <- unique(network[["target"]])

    targets <- rownames(mat)
    conditions <- colnames(mat)

    # Extract matrix of weights
    weight_mat <- network %>%
        pivot_wider_profile(
            id_cols = .data$source,
            names_from = .data$target,
            values_from = .data$weight,
            to_matrix = TRUE,
            to_sparse = sparse,
            values_fill = 0
        )

    weight_mat <- as.matrix(weight_mat)

    # This fixes the wrong denominator defined in contribution
    weight_mat <- weight_mat/rowSums(abs(weight_mat))

    # Analysis ----------------------------------------------------------------
    withr::with_seed(seed, {
        .wmean_analysis(mat, weight_mat, shared_targets, times, randomize_type)
    })

}

# Helper functions --------------------------------------------------------

#' Wrapper to execute run_wmean() logic once finished preprocessing of data
#'
#' @inherit run_wmean description
#'
#' @inheritParams run_wmean
#' @param weight_mat Matrix that corresponds to the multiplication of the mor
#'  column with likelihood divided over the contribution.
#' @param shared_targets Target nodes that are shared between the
#'  `mat` and `network`.
#'
#' @inherit run_wmean return
#'
#' @keywords internal
#' @noRd
.wmean_analysis <- function(mat, weight_mat, shared_targets, times, randomize_type) {
    # Thus, it is only necessary to define if we want
    # to evaluate a random model or not.
    wmean_run <- partial(
        .wmean_run,
        mat = mat,
        weight_mat = weight_mat,
        shared_targets = shared_targets,
        randomize_type = randomize_type
    )

    # Run model for random data
    map_dfr(seq_len(times), ~ wmean_run(random = TRUE)) %>%
        group_by(.data$source, .data$condition) %>%
        summarise(
            null_distribution = list(.data$value),
            null_mean = mean(.data$value),
            null_sd = stats::sd(.data$value),
            .groups = "drop"
        ) %>%
        # Run the true model and joined to random.
        left_join(y = wmean_run(random = FALSE), by = c("source", "condition")) %>%
        # Calculate scores
        mutate(
            z_score = (.data$value - .data$null_mean) / .data$null_sd,
            z_score = replace_na(.data$z_score, 0),
            p_value = map2_dbl(
                .x = .data$null_distribution,
                .y = .data$value,
                .f = ~ sum(abs(.x) > abs(.y)) / length(.x)
            ),
            # Limit empirical p-value to lower bound 1/times and upper bound
            # (times-1)/times
            p_value = if_else(.data$p_value == 0, 1/times, .data$p_value),
            p_value = if_else(.data$p_value == 1, (times-1)/times, .data$p_value),
            p_value = if_else(.data$p_value >= 0.5, 1-.data$p_value, .data$p_value),
            p_value = .data$p_value * 2,
            c_score = .data$value * (-log10(.data$p_value))
        ) %>%
        # Reformat results
        select(-contains("null")) %>%
        rename(corr_wmean = .data$c_score, wmean = .data$value, norm_wmean = .data$z_score) %>%
        pivot_longer(
            cols = c(.data$corr_wmean, .data$wmean, .data$norm_wmean),
            names_to = "statistic",
            values_to = "score"
        ) %>%
        arrange(.data$statistic, .data$source, .data$condition) %>%
        select(.data$statistic, .data$source, .data$condition, .data$score, .data$p_value)
}

#' Wrapper to run wmean one time
#'
#' @inheritParams .wmean_analysis
#' @inherit .wmean_evaluate_model return
#' @keywords internal
#' @noRd
.wmean_run <- function(mat, weight_mat, shared_targets, random, randomize_type) {
    .wmean_map_model_data(mat, shared_targets, random, randomize_type) %>%
        .wmean_evaluate_model(weight_mat)
}

#' Calculate mean weight
#'
#' @inheritParams .wmean_analysis
#' @keywords internal
#' @noRd
.wmean_calculate_weight <- function(network) {
    network %>%
        add_count(.data$source, name = "contribution") %>%
        transmute(
            .data$source,
            .data$target,
            weight = .data$mor * .data$likelihood
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
#' @inheritParams .wmean_analysis
#' @keywords internal
#' @noRd
.wmean_map_model_data <- function(mat, shared_targets, random, randomize_type) {
    if (random) {
        randomize_matrix(mat, randomize_type = randomize_type)[shared_targets, , drop = FALSE]
    } else {
        mat[shared_targets, , drop = FALSE]
    }
}

#' Evaluate model
#'
#' The evaluation model consists of evaluating the multiplication of the
#' weights by the factor of interest and comparing it against results
#' from permutations of the matrix of values of interest.
#'
#' @inheritParams .wmean_analysis
#'
#' @return A dataframe with three columns:
#'  source (source nodes), condition (colnames of mat) and value (score).
#'
#' @keywords internal
#' @noRd
.wmean_evaluate_model <- function(mat, weight_mat) {
    (weight_mat %*% mat) %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column("source") %>%
        pivot_longer(-.data$source, names_to = "condition")
}
