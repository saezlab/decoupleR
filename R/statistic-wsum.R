#' Weighted Sum (WSUM)
#'
#' @description
#' Calculates regulatory activities using WSUM.
#'
#' @details
#' WSUM infers regulator activities by first multiplying each target feature by
#' its associated weight which then are summed to an enrichment score
#' `wsum`. Furthermore, permutations of random target features can be
#' performed to obtain a null distribution that can be used to compute a z-score
#' `norm_wsum`, or a corrected estimate `corr_wsum` by multiplying
#' `wsum` by the minus log10 of the obtained empirical p-value.
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
#' run_wsum(mat, net, minsize=0)
run_wsum <- function(mat,
                     network,
                     .source = source,
                     .target = target,
                     .mor = mor,
                     .likelihood = likelihood,
                     times = 100,
                     seed = 42,
                     sparse = TRUE,
                     randomize_type = "rows",
                     minsize = 5
                     ) {

    # NSE vs. R CMD check workaround
    c_score <- condition <- corr_wsum <- likelihood <- mor <- norm_wsum <-    null_distribution <- null_mean <- null_sd <- p_value <- score <- source <-    statistic <- target <- value <- weight <- wsum <- z_score <- NULL

    # Before to start ---------------------------------------------------------
    if (times < 2) {
        rlang::abort(message = stringr::str_glue("Parameter 'times' must be greater than or equal to 2, but {times} was passed."))
    }

    # Check for NAs/Infs in mat
    mat <- check_nas_infs(mat)

    network <- network %>%
        rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
    network <- filt_minsize(rownames(mat), network, minsize)

    # Preprocessing -----------------------------------------------------------

    # Calculate the weights that will be used for the evaluation of the model
    network <- network %>%
        filter(target %in% rownames(mat)) %>%
        .wsum_calculate_weight()

    # Extract labels that will map to the expression and profile matrices
    shared_targets <- unique(network[["target"]])

    targets <- rownames(mat)
    conditions <- colnames(mat)

    # Extract matrix of weights
    weight_mat <- network %>%
        pivot_wider_profile(
            id_cols = source,
            names_from = target,
            values_from = weight,
            to_matrix = TRUE,
            to_sparse = sparse,
            values_fill = 0
        )

    weight_mat <- as.matrix(weight_mat)

    # Analysis ----------------------------------------------------------------
    withr::with_seed(seed, {
        .wsum_analysis(mat, weight_mat, shared_targets, times, randomize_type)
    })

}

# Helper functions --------------------------------------------------------

#' Wrapper to execute run_wsum() logic once finished preprocessing of data
#'
#' @inherit run_wsum description
#'
#' @inheritParams run_wsum
#' @param weight_mat Matrix that corresponds to the multiplication of the mor
#'  column with likelihood divided over the contribution.
#' @param shared_targets Target nodes that are shared between the
#'  `mat` and `network`.
#'
#' @inherit run_wsum return
#'
#' @keywords internal
#' @noRd
.wsum_analysis <- function(mat, weight_mat, shared_targets, times, randomize_type) {
    # Thus, it is only necessary to define if we want
    # to evaluate a random model or not.
    wsum_run <- partial(
        .wsum_run,
        mat = mat,
        weight_mat = weight_mat,
        shared_targets = shared_targets,
        randomize_type = randomize_type
    )


    # NSE vs. R CMD check workaround
    corr_wsum <- wsum <- value <- source <- condition <- p_value <-
    null_mean <- null_distribution <- null_sd <- z_score <- norm_wsum <-
    statistic <- score <- NULL

    # Run model for random data
    map_dfr(seq_len(times), ~ wsum_run(random = TRUE)) %>%
        group_by(source, condition) %>%
        summarise(
            null_distribution = list(value),
            null_mean = mean(value),
            null_sd = stats::sd(value),
            .groups = "drop"
        ) %>%
        # Run the true model and joined to random.
        left_join(y = wsum_run(random = FALSE), by = c("source", "condition")) %>%
        # Calculate scores
        mutate(
            z_score = (value - null_mean) / null_sd,
            z_score = replace_na(z_score, 0),
            p_value = map2_dbl(
                .x = null_distribution,
                .y = value,
                .f = ~ sum(abs(.x) > abs(.y)) / length(.x)
            ),
            # Limit empirical p-value to lower bound 1/times and upper bound
            # (times-1)/times
            p_value = if_else(p_value == 0, 1/times, p_value),
            p_value = if_else(p_value == 1, (times-1)/times, p_value),
            p_value = if_else(p_value >= 0.5, 1-p_value, p_value),
            p_value = p_value * 2,
            c_score = value * (-log10(p_value))
        ) %>%
        # Reformat results
        select(-contains("null")) %>%
        rename(corr_wsum = c_score, wsum = value, norm_wsum = z_score) %>%
        pivot_longer(
            cols = c(corr_wsum, wsum, norm_wsum),
            names_to = "statistic",
            values_to = "score"
        ) %>%
        arrange(statistic, source, condition) %>%
        select(statistic, source, condition, score, p_value)
}

#' Wrapper to run wsum one time
#'
#' @inheritParams .wsum_analysis
#' @inherit .wsum_evaluate_model return
#' @keywords internal
#' @noRd
.wsum_run <- function(mat, weight_mat, shared_targets, random, randomize_type) {
    .wsum_map_model_data(mat, shared_targets, random, randomize_type) %>%
        .wsum_evaluate_model(weight_mat)
}

#' Calculate sum weight
#'
#' @inheritParams .wsum_analysis
#' @keywords internal
#' @noRd
.wsum_calculate_weight <- function(network) {
    network %>%
        add_count(source, name = "contribution") %>%
        transmute(
            source,
            target,
            weight = mor * likelihood
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
#' @inheritParams .wsum_analysis
#' @keywords internal
#' @noRd
.wsum_map_model_data <- function(mat, shared_targets, random, randomize_type) {
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
#' @inheritParams .wsum_analysis
#'
#' @return A dataframe with three columns:
#'  source (source nodes), condition (colnames of mat) and value (score).
#'
#' @keywords internal
#' @noRd
.wsum_evaluate_model <- function(mat, weight_mat) {
    (weight_mat %*% mat) %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column("source") %>%
        pivot_longer(-source, names_to = "condition")
}
