#' Multivariate Linear Model (MLM)
#'
#' @description
#' Calculates regulatory activities by fitting multivariate linear models (MLM)
#'
#' @details
#' MLM fits a multivariate linear model to estimate regulatory activities.
#' MLM transforms a given network into an adjacency matrix, placing sources as
#' columns and targets as rows. The matrix is filled with the associated weights
#' for each interaction. This matrix is used to fit a linear model to predict
#' the observed molecular readouts per sample. The obtained t-values from the
#' fitted model are the activities of the regulators.
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
#' run_mlm(mat, network, .source='tf')
run_mlm <- function(mat,
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
    convert_to_mlm({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Preprocessing -----------------------------------------------------------
  .mlm_preprocessing(network, mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
  {
    .mlm_analysis(.$mat, .$mor_mat)
  }
}

# Helper functions ------------------------------------------------------
#' mlm preprocessing
#'
#' - Get only the intersection of target genes between `mat` and `network`.
#' - Transform tidy `network` into `matrix` representation with `mor` as value.
#' - If `center` is true, then the expression values are centered by the
#'   mean of expression across the conditions.
#'
#' @inheritParams run_mlm
#'
#' @return A named list of matrices to evaluate in `.mlm_analysis()`.
#'  - mat: Genes as rows and conditions as columns.
#'  - mor_mat: Genes as rows and columns as source.
#' @keywords intern
#' @noRd
.mlm_preprocessing <- function(network, mat, center, na.rm, sparse) {
  shared_targets <- intersect(
    rownames(mat),
    network$target
  )

  mat <- mat[shared_targets, , drop=FALSE]

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
    .[shared_targets, , drop=FALSE]

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
    .[shared_targets, , drop=FALSE]

  weight_mat <- mor_mat * likelihood_mat

  if (center) {
    mat <- mat - rowMeans(mat, na.rm)
  }

  list(mat = mat, mor_mat = weight_mat)
}

#' Wrapper to execute run_mlm() logic one finished preprocessing of data
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @inheritParams run_mlm
#' @param mor_mat
#'
#' @inherit run_mlm return
#' @keywords intern
#' @noRd
.mlm_analysis <- function(mat, mor_mat) {
  mlm_evaluate_model <- partial(
    .f = .mlm_evaluate_model,
    mat = mat,
    mor_mat = mor_mat
  )

  # Allocate the space for all conditions and evaluate the proposed model.
  expand_grid(
    condition = colnames(mat)
  ) %>%
    rowwise(.data$condition) %>%
    mutate(model = list(mlm_evaluate_model(.data$condition)), statistic='mlm') %>%
    unnest(.data$model) %>%
    select(.data$statistic, .data$source, .data$condition, .data$score, .data$p_value)
}

#' Wrapper to run mlm per sample (condition) at time
#'
#' @keywords internal
#' @noRd
.mlm_evaluate_model <- function(condition, mat, mor_mat) {
  #data <- cbind(data.frame(y=mat[ , condition]), mor_mat)
  fit <- lm(mat[ , condition] ~ mor_mat) %>%
      summary()
  scores <- as.vector(fit$coefficients[,3][-1])
  pvals <- as.vector(fit$coefficients[,4][-1])
  tibble(score=scores, p_value=pvals, source=colnames(mor_mat))
}
