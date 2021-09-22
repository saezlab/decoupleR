#' Univariate Linear Model (ULM)
#'
#' @description
#' Calculates regulatory activities by fitting univariate linear models (ULM)
#'
#' @details
#' ULM fits a (univariate) linear model to estimate regulatory activities. ULM
#' fits a linear model that predicts the observed molecular using the given
#' weights of a regulon as a single co-variate. The obtained t-value
#' from the fitted model is the activity of the regulon. THis approach was first
#' described in:
#' [Improved detection of tumor suppressor events in single-cell RNA-Seq data](
#' https://www.nature.com/articles/s41525-020-00151-y?elqTrackId=d7efb03cf5174fe2ba84e1c34d602b13).
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Logical value indicating if the generated profile matrix
#'  should be sparse.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
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
#' run_ulm(mat, network, .source='tf')
run_ulm <- function(mat,
                      network,
                      .source = .data$source,
                      .target = .data$target,
                      .mor = .data$mor,
                      .likelihood = .data$likelihood,
                      sparse = FALSE,
                      center = FALSE,
                      na.rm = FALSE) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    # Convert to standard tibble: source-target-mor.
    network <- network %>%
        convert_to_ulm({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

    # Preprocessing -----------------------------------------------------------
    .ulm_preprocessing(network, mat, center, na.rm, sparse) %>%
        # Model evaluation --------------------------------------------------------
        {
            .ulm_analysis(.$mat, .$mor_mat)
        }
}

# Helper functions ------------------------------------------------------
#' ulm preprocessing
#'
#' - Get only the intersection of target genes between `mat` and `network`.
#' - Transform tidy `network` into `matrix` representation with `mor` as value.
#' - If `center` is true, then the expression values are centered by the
#'   mean of expression across the conditions.
#'
#' @inheritParams run_ulm
#'
#' @return A named list of matrices to evaluate in `.ulm_analysis()`.
#'  - mat: Genes as rows and conditions as columns.
#'  - mor_mat: Genes as rows and columns as source.
#' @keywords intern
#' @noRd
.ulm_preprocessing <- function(network, mat, center, na.rm, sparse) {
    shared_targets <- intersect(
        rownames(mat),
        network$target
    )

    mat <- mat[shared_targets, ]

    mor_mat <- network %>%
        filter(.data$target %in% shared_targets) %>%
        pivot_wider_profile(
            id_cols = .data$target,
            names_from = .data$source,
            values_from = .data$mor,
            values_fill = 0,
            to_matrix = TRUE,
            to_sparse = sparse
        ) %>%
        .[shared_targets, ]

    likelihood_mat <- network %>%
        filter(.data$target %in% shared_targets) %>%
        pivot_wider_profile(
            id_cols = .data$target,
            names_from = .data$source,
            values_from = .data$likelihood,
            values_fill = 0,
            to_matrix = TRUE,
            to_sparse = sparse
        ) %>%
        .[shared_targets, ]

    weight_mat <- mor_mat * likelihood_mat

    if (center) {
        mat <- mat - rowMeans(mat, na.rm)
    }

    list(mat = mat, mor_mat = weight_mat)
}

#' Wrapper to execute run_ulm() logic one finished preprocessing of data
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @inheritParams run_ulm
#' @param mor_mat
#'
#' @inherit run_ulm return
#' @keywords intern
#' @noRd
.ulm_analysis <- function(mat, mor_mat) {
    ulm_evaluate_model <- partial(
        .f = .ulm_evaluate_model,
        mat = mat,
        mor_mat = mor_mat
    )

    # Allocate the space for all combinations of sources and conditions
    # and evaluate the proposed model.
    expand_grid(
        source = colnames(mor_mat),
        condition = colnames(mat)
    ) %>%
        rowwise(.data$source, .data$condition) %>%
        summarise(
            score = ulm_evaluate_model(.data$source, .data$condition),
            .groups = "drop"
        ) %>%
        transmute(statistic = "ulm", .data$source, .data$condition, .data$score)
}

#' Wrapper to run ulm one source (source) per sample (condition) at time
#'
#' @keywords internal
#' @noRd
.ulm_evaluate_model <- function(source, condition, mat, mor_mat) {
        speedlm.fit(
            y = mat[, condition],
            X = cbind(1, mor_mat[, source])
        ) %>%
            summary() %>%
            pluck("coefficients", "t", 2, .default = NA)
}
