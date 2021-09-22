#' Univariate Decision Tree (UDT)
#'
#' UDT fits a (univariate) decision tree to estimate regulatory activities. UDT
#' fits a decision tree that predicts the observed molecular using the given
#' weights of a regulon as a single co-variate. The obtained feature importance
#' from the fitted model is the activity of the regulon.
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
#' @import rpart
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_udt(mat, network, .source='tf')
run_udt <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    .mor = .data$mor,
                    .likelihood = .data$likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    min_n = 2,
                    seed = 42
) {
  set.seed(seed)
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  # Before to start ---------------------------------------------------------
  # Convert to standard tibble: source-target-mor.
  network <- network %>%
    convert_to_ulm({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Preprocessing -----------------------------------------------------------
  .udt_preprocessing(network, mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
  {
    .udt_analysis(.$mat, .$mor_mat, min_n, seed)
  }
}

# Helper functions ------------------------------------------------------
#' udt preprocessing
#'
#' - Get only the intersection of target genes between `mat` and `network`.
#' - Transform tidy `network` into `matrix` representation with `mor` as value.
#' - If `center` is true, then the expression values are centered by the
#'   mean of expression across the conditions.
#'
#' @inheritParams run_udt
#'
#' @return A named list of matrices to evaluate in `.udt_analysis()`.
#'  - mat: Genes as rows and conditions as columns.
#'  - mor_mat: Genes as rows and columns as source.
#' @keywords intern
#' @noRd
.udt_preprocessing <- function(network, mat, center, na.rm, sparse) {
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

#' Wrapper to execute run_udt() logic once finished preprocessing of data
#'
#'
#' @inheritParams run_udt
#' @param mor_mat
#'
#' @inherit run_udt return
#' @keywords intern
#' @noRd
.udt_analysis <- function(mat, mor_mat, min_n, seed) {
  udt_evaluate_model <- partial(
    .f = .udt_evaluate_model,
    mat = mat,
    mor_mat = mor_mat,
    min_n = min_n
  )

  # Allocate the space for all conditions and evaluate the proposed model.
  temp <- expand_grid(
    source = colnames(mor_mat),
    condition = colnames(mat)
  )

  score <- seq_len(nrow(temp)) %>%
    map_dbl(~udt_evaluate_model(temp %>% pluck("source", .x),
                                temp %>% pluck("condition", .x)))

  bind_cols(temp, score = score) %>%
    transmute(statistic = "udt", .data$source, .data$condition, .data$score)
}

#' Wrapper to run udt per a sample (condition) at time
#'
#' @keywords internal
#' @noRd
.udt_evaluate_model <- function(source, condition, mat, mor_mat, min_n) {
  data <- tibble(x = mat[, condition, drop=F] , y = mor_mat[, source])
  score <- rpart::rpart(y~x, data, minsplit=min_n) %>% pluck("variable.importance")

  if (is.null(score)) {
    score <- 0
    names(score) <- source
  }
  score
}
