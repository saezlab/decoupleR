#' PSCIRA (Permutation Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculate the regulatory activity of each source by multiplying the expression
#' values of its objectives with their corresponding associated profiles for
#' each given condition.The result is equal to the z-score of the found value
#' compared to its null distribution.
#'
#' @details
#' PSCIRA estimates the regulatory activity by performing a linear combination
#' of the expression profile and the target profile of a given source. In the
#' target profile of a source, any regulon member is assigned a +1 for activating
#' interactions and a -1 for inhibitory interactions. Additionally, any
#' positive weight can be added to each interaction. A normalized score per
#' score is calculated by calculating a z-score using a null distribution
#' computed from n permutations (times parameter). A corrected score per score is
#' calculated by calculating a p-value using a null distribution computed from n
#' permutations.
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
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_pscira(mat, network, .source='tf')
run_pscira <- function(mat,
                       network,
                       .source = .data$source,
                       .target = .data$target,
                       .mor = .data$mor,
                       .likelihood = .data$likelihood,
                       sparse = TRUE,
                       times = 100,
                       seed = 42) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    if (times < 2) {
        rlang::abort(message = stringr::str_glue("Parameter 'times' must be greater than or equal to 2, but {times} was passed."))
    }

    # Preprocessing -----------------------------------------------------------

    # Convert to standard tibble: source-target-mor.
    network <- network %>%
        convert_to_pscira({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

    # Extract labels that will map to the expression and profile matrices
    sources <- network %>%
        pull(.data$source) %>%
        unique()

    # Ensures column matching, expands the target profile to encompass
    # all targets in the expression matrix for each source, and converts
    # the result to a matrix.
    mor_mat <- network %>%
        get_profile_of(
            sources = list(source = sources, target = rownames(mat)),
            values_fill = list(mor = 0)
        ) %>%
        pivot_wider_profile(
            id_cols = .data$source,
            names_from = .data$target,
            values_from = .data$mor,
            to_matrix = TRUE,
            to_sparse = sparse
        )

    likelihood_mat <- network %>%
        get_profile_of(
            sources = list(source = sources, target = rownames(mat)),
            values_fill = list(likelihood = 0)
        ) %>%
        pivot_wider_profile(
            id_cols = .data$source,
            names_from = .data$target,
            values_from = .data$likelihood,
            to_matrix = TRUE,
            to_sparse = sparse
        )

    weight_mat <- mor_mat * likelihood_mat

    # Convert to matrix to ensure that matrix multiplication works
    # in case mat is a labeled dataframe.
    mat <- as.matrix(mat)

    # Evaluate model ----------------------------------------------------------
    withr::with_seed(seed, {
        .pscira_analysis(mat, weight_mat, times)
    })
}

# Helper functions --------------------------------------------------------

#' Wrapper to execute run_pscira() logic one finished preprocessing of data
#'
#' @inheritParams run_pscira
#' @param mor_mat Matrix that corresponds to the mor of the
#' target genes (columns) of a source (rows).
#'
#' @inherit run_pscira return
#' @keywords intern
#' @noRd
.pscira_analysis <- function(mat, mor_mat, times) {
    pscira_run <- partial(
        .f = .pscira_run,
        mat = mat,
        mor_mat = mor_mat
    )

    map_dfr(seq_len(times), ~ pscira_run(random = TRUE)) %>%
        group_by(.data$source, .data$condition) %>%
        summarise(
            null_distribution = list(.data$value),
            .mean = mean(.data$value),
            .sd = sd(.data$value),
            .groups = "drop"
        ) %>%
        left_join(pscira_run(random = FALSE), by = c("source", "condition")) %>%
        mutate(
            score = (.data$value - .data$.mean) / .data$.sd,
            score = replace_na(.data$score, 0),
            p_value = map2_dbl(
                .x = .data$null_distribution,
                .y = .data$value,
                .f = ~ sum(abs(.x) > abs(.y)) / length(.x)
            ),
            c_p_value = ifelse(.data$p_value == 0, 1/length(.data$null_distribution), .data$p_value),
            c_score = .data$value * (-log10(.data$c_p_value))
        ) %>%
        rename(corrected_pscira = .data$c_score, normalized_pscira = .data$score, pscira = .data$value) %>%
        pivot_longer(
            cols = c(.data$corrected_pscira, .data$normalized_pscira, .data$pscira),
            names_to = "statistic",
            values_to = "score"
        ) %>%
        arrange(.data$statistic, .data$source, .data$condition) %>%
        select(.data$statistic, .data$source, .data$condition, .data$score, .data$p_value)
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
#' Calculates the regulatory activity of all sources with respect to its
#' associated profile for each condition.
#'
#' @inheritParams .pscira_run
#'
#' @return Tibble with source activity for each source-sample pair.
#' @keywords internal
#' @noRd
.pscira_evaluate_model <- function(mat, mor_mat) {
    (mor_mat %*% mat) %>%
        as.matrix() %>%
        as.data.frame() %>%
        rownames_to_column("source") %>%
        pivot_longer(-.data$source, names_to = "condition")
}
