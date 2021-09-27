#' Multivariate Decision Trees (MDT)
#'
#' @description
#' Calculates regulatory activities by fitting multivariate desicion trees (MDT)
#'
#' @details
#' MDT fits a multivariate ensemble of decision trees (random forest) to
#' estimate regulatory activities. MDT transforms a given network into an
#' adjacency matrix, placing sources as columns and targets as rows. The matrix
#' is filled with the associated weights for each interaction. This matrix is
#' used to fit a random forest model to predict the observed molecular
#' readouts per sample. The obtained feature importances from the fitted model
#' are the activities of the regulators.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Logical value indicating if the generated profile matrix
#'  should be sparse.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param trees An integer for the number of trees contained in the ensemble.
#' @param min_n An integer for the minimum number of data points in a node that
#' are required for the node to be split further.
#' @param nproc Number of cores to use for computation.
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
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
#' @import ranger
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_mdt(mat, network, .source='tf')
run_mdt <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    .mor = .data$mor,
                    .likelihood = .data$likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    trees = 10,
                    min_n = 20,
                    nproc = 4,
                    seed = 42
) {
  set.seed(seed)
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  # Before to start ---------------------------------------------------------
  # Convert to standard tibble: source-target-mor.
  network <- network %>%
    convert_to_mlm({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Preprocessing -----------------------------------------------------------
  .mdt_preprocessing(network, mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
  {
    .mdt_analysis(.$mat, .$mor_mat, trees, min_n, nproc)
  }
}

# Helper functions ------------------------------------------------------
#' mdt preprocessing
#'
#' - Get only the intersection of target genes between `mat` and `network`.
#' - Transform tidy `network` into `matrix` representation with `mor` as value.
#' - If `center` is true, then the expression values are centered by the
#'   mean of expression across the conditions.
#'
#' @inheritParams run_mdt
#'
#' @return A named list of matrices to evaluate in `.nolea_analysis()`.
#'  - mat: Genes as rows and conditions as columns.
#'  - mor_mat: Genes as rows and columns as source.
#' @keywords intern
#' @noRd
.mdt_preprocessing <- function(network, mat, center, na.rm, sparse) {
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

#' Wrapper to execute run_mdt() logic one finished preprocessing of data
#'
#'
#' @inheritParams run_mdt
#' @param mor_mat
#'
#' @inherit run_mdt return
#' @keywords intern
#' @noRd
.mdt_analysis <- function(mat, mor_mat, trees, min_n, nproc) {
  mdt_evaluate_model <- partial(
    .f = .mdt_evaluate_model,
    mat = mat,
    mor_mat = mor_mat,
    trees = trees,
    min_n = min_n,
    nproc = nproc
  )

  # Allocate the space for all conditions and evaluate the proposed model.
  expand_grid(
    condition = colnames(mat)
  ) %>%
    rowwise(.data$condition) %>%
    summarise(
      score = mdt_evaluate_model(.data$condition),
      source = colnames(mor_mat),
      .groups = "drop"
    ) %>%
    transmute(statistic = "mdt", .data$source, .data$condition, .data$score
    ) %>%
    arrange(source)
}

#' Wrapper to run mdt per a sample (condition) at time
#'
#' @keywords internal
#' @noRd
.mdt_evaluate_model <- function(condition, mat, mor_mat, trees, min_n, nproc) {
  ranger::ranger(condition ~ ., data = data.frame(condition=mat[, condition], mor_mat),
                 num.trees = trees,
                 importance = "impurity",
                 min.node.size = min_n,
                 num.threads = nproc) %>%
    pluck("variable.importance")

}
