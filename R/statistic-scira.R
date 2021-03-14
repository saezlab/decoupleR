#' SCIRA (Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculates TF activity according to
#' [Improved detection of tumor suppressor events in single-cell RNA-Seq data](
#' https://www.nature.com/articles/s41525-020-00151-y?elqTrackId=d7efb03cf5174fe2ba84e1c34d602b13)
#' .
#'
#' @details
#' Estimation of regulatory activity: A linear regression of the expression
#' profile is performed against the "target profile" of the given TF, where
#' in the target profile, any regulon member is assigned a `+1` for activating
#' interactions and a `-1` for inhibitory interactions. All other genes not
#' members of the TF's regulon are assigned a value o `0`. TF activity is then
#' defined as the t-statistic of this linear regression.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Logical value indicating if the generated profile matrix
#'  should be sparse.
#' @param fast Logical value indicating if the lineal model must be calculated
#' with [speedglm::speedlm.fit()] or with base [stats::lm()].
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats coef lm summary.lm
#' @importFrom speedglm speedlm.fit
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_scira(mat, network, tf, target, mor)
run_scira <- function(
    mat,
    network,
    .source = .data$tf,
    .target = .data$target,
    .mor = .data$mor,
    sparse = FALSE,
    fast = TRUE,
    center = TRUE,
    na.rm = FALSE) {

    # Before to start ---------------------------------------------------------
    # Convert to standard tibble: tf-target-mor.
    network <- network %>%
        convert_to_scira({{ .source }}, {{ .target }}, {{ .mor }})

    # Preprocessing -----------------------------------------------------------
    .scira_preprocessing(network, mat, center, na.rm, sparse) %>%

    # Model evaluation --------------------------------------------------------
    {
        .scira_analysis(.$mat, .$mor_mat, fast)
    }
}

# Helper functions ------------------------------------------------------
#' Scira preprocessing
#'
#' - Get only the intersection of target genes between `mat` and `network`.
#' - Transform tidy `network` into `matrix` representation with `mor` as value.
#' - If `center` is true, then the expression values are centered by the
#'   mean of expression across the conditions.
#'
#' @inheritParams run_scira
#'
#' @return A named list of matrices to evaluate in `.scira_analysis()`.
#'  - mat: Genes as rows and conditions as columns.
#'  - mor_mat: Genes as rows and columns as tfs.
#' @keywords intern
#' @noRd
.scira_preprocessing <- function(network, mat, center, na.rm, sparse) {

    shared_targets <- intersect(
        rownames(mat),
        network$target
    )

    mat <- mat[shared_targets, ]

    mor_mat <- network %>%
        filter(.data$target %in% shared_targets) %>%
        pivot_wider_profile(
            id_cols = .data$target,
            names_from = .data$tf,
            values_from = .data$mor,
            values_fill = 0,
            to_matrix = TRUE,
            to_sparse = sparse
        ) %>%
        .[shared_targets, ]

    if (center) {
        mat <- mat - rowMeans(mat, na.rm)
    }

    list(mat = mat, mor_mat = mor_mat)
}

#' Wrapper to execute run_scira() logic one finished preprocessing of data
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @inheritParams run_scira
#' @param mor_mat
#'
#' @inherit run_scira return
#' @keywords intern
#' @noRd
.scira_analysis <- function(mat, mor_mat, fast) {
    scira_evaluate_model <- partial(
        .f = .scira_evaluate_model,
        mat = mat,
        mor_mat = mor_mat,
        fast = fast
    )

    # Allocate the space for all combinations of sources and conditions
    # and evaluate the proposed model.
    expand_grid(
        tf = colnames(mor_mat),
        condition = colnames(mat)
    ) %>%
        rowwise(.data$tf, .data$condition) %>%
        summarise(
            score = scira_evaluate_model(.data$tf, .data$condition),
            .groups = "drop"
        ) %>%
        transmute(statistic = "scira", .data$tf, .data$condition, .data$score)
}

#' Wrapper to run scira one tf (source) per sample (condition) at time
#'
#' @keywords internal
#' @noRd
.scira_evaluate_model <- function(source, condition, mat, mor_mat, fast) {
    if (fast) {
        speedlm.fit(
            y = mat[, condition],
            X = cbind(1, mor_mat[, source])
        ) %>%
            summary() %>%
            pluck("coefficients", "t", 2, .default = NA)
    } else {
        lm(mat[, condition] ~ mor_mat[, source]) %>%
            summary() %>%
            coefficients() %>%
            .[, "t value"] %>%
            pluck(2, .default = NA)
    }
}
