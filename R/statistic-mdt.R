#' MDT (Multivariate Decision Tree)
#'
#' @description
#'
#'
#' @details
#'
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
                    num.threads = 4,
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
    .mdt_analysis(.$mat, .$mor_mat, trees, num.threads)
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
.mdt_analysis <- function(mat, mor_mat, trees, num.threads) {
  mdt_evaluate_model <- partial(
    .f = .mdt_evaluate_model,
    mat = mat,
    mor_mat = mor_mat,
    trees = trees,
    num.threads = num.threads
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
.mdt_evaluate_model <- function(condition, mat, mor_mat, trees, num.threads) {
  ranger::ranger(condition ~ ., data = data.frame(condition=mat[, condition], mor_mat), importance = "impurity", num.threads = num.threads) %>% 
    pluck("variable.importance")
  
}