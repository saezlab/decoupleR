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
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_pscira(mat, network, tf, target, mor)
run_pscira <- function(mat,
                       network,
                       .source = .data$tf,
                       .target = .data$target,
                       .mor = .data$mor,
                       .likelihood = .data$likelihood,
                       sparse = TRUE,
                       times = 10,
                       seed = 42) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    if (times < 2) {
        rlang::abort(message = stringr::str_glue("Parameter 'times' must be greater than or equal to 2, but {times} was passed."))
    }

    # Preprocessing -----------------------------------------------------------

    # Convert to standard tibble: tf-target-mor.
    network <- network %>%
        convert_to_pscira({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

    # Extract labels that will map to the expression and profile matrices
    tfs <- network %>%
        pull(.data$tf) %>%
        unique()

    # Ensures column matching, expands the target profile to encompass
    # all targets in the expression matrix for each source, and converts
    # the result to a matrix.
    mor_mat <- network %>%
        get_profile_of(
            sources = list(tf = tfs, target = rownames(mat)),
            values_fill = list(mor = 0)
        ) %>%
        pivot_wider_profile(
            id_cols = .data$tf,
            names_from = .data$target,
            values_from = .data$mor,
            to_matrix = TRUE,
            to_sparse = sparse
        )

    likelihood_mat <- network %>%
        get_profile_of(
            sources = list(tf = tfs, target = rownames(mat)),
            values_fill = list(likelihood = 0)
        ) %>%
        pivot_wider_profile(
            id_cols = .data$tf,
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
#' target genes (columns) of a tf (rows).
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
        group_by(.data$tf, .data$condition) %>%
        summarise(
            .mean = mean(.data$value),
            .sd = sd(.data$value),
            .groups = "drop"
        ) %>%
        left_join(pscira_run(random = FALSE), by = c("tf", "condition")) %>%
        mutate(
            score = (.data$value - .data$.mean) / .data$.sd,
            score = replace_na(.data$score, 0)
        ) %>%
        rename(normalized_pscira = .data$score, pscira = .data$value) %>%
        pivot_longer(
            cols = c(.data$normalized_pscira, .data$pscira),
            names_to = "statistic",
            values_to = "score"
        ) %>%
        arrange(.data$statistic, .data$tf, .data$condition) %>%
        select(.data$statistic, .data$tf, .data$condition, .data$score)
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
