#' Multivariate Linear Model (MLM) with regularization
#'
#' @description
#' Calculates regulatory activities by fitting multivariate linear models (MLM)
#'
#'TODO: will modify the manual
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
#' @param sparse Deprecated parameter.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param minsize Integer indicating the minimum number of targets per source.
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
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_mlm(mat, network, .source='tf')
run_mlmreg <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    .mor = .data$mor,
                    .likelihood = .data$likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    minsize = 5,
                    alpha = 1,
                    low.limits = -Inf) {

  # "norm" parameter can only be L1 or L2
  stopifnot((alpha >= 0 ) && (alpha <= 1))

  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  # Before to start ---------------------------------------------------------
  # Convert to standard tibble: source-target-mor.
  network <- network %>%
    rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
  network <- filt_minsize(rownames(mat), network, minsize)

  # Preprocessing -----------------------------------------------------------
  .fit_preprocessing(network, mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
  {
    .mlmreg_analysis(.$mat, .$mor_mat, alpha)
  } %>%
    ungroup()
}

#' Wrapper to execute run_mlmreg() logic one finished preprocessing of data
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @inheritParams run_mlmreg
#' @param mor_mat
#'
#' @inherit run_mlmreg return
#' @keywords intern
#' @noRd
.mlmreg_analysis <- function(mat, mor_mat, alpha) {
  mlmreg_evaluate_model <- partial(
    .f = .mlmreg_evaluate_model,
    mat = mat,
    mor_mat = mor_mat,
    alpha = alpha
  )

  # Allocate the space for all conditions and evaluate the proposed model.
  expand_grid(
    condition = colnames(mat)
  ) %>%
    rowwise(.data$condition) %>%
    mutate(model = list(mlmreg_evaluate_model(.data$condition)), statistic='mlmreg') %>%
    unnest(.data$model) %>%
    select(.data$statistic, .data$source, .data$condition, .data$score)
}

#' Wrapper to run mlmreg per sample (condition) at time
#'
#' @keywords internal
#' @noRd
.mlmreg_evaluate_model <- function(condition, mat, mor_mat, alpha=1) {

  fit <- glmnet(mor_mat, mat[, condition],
                lambda.min.ratio=0.0001, nlambda=100,
                lower.limits = -Inf,
                alpha=alpha, # control ridge vs lasso
                standardize=F, family='gaussian')

  coef = as.matrix(fit$beta)
  coef = round(coef, digits=6)
  colnames(coef) = round(fit$lambda, digits = 9)

  subcol = abs(colSums(coef)) >0
  coef1 = coef[, subcol]
  # norm_beta = scale(abs(coef1)**(2-alpha), center=FALSE,
                    # scale=colSums(abs(coef1)**(2-alpha))/fit$dev.ratio[subcol] ) # adjust statistics by R^2, not helping performance
  norm_beta = scale(abs(coef1)**(2-alpha), center=FALSE,
                    scale=colSums(abs(coef1)**(2-alpha)))
  beta_max = apply(norm_beta, 1, max)

  beta_max_idx = apply(norm_beta, 1, which.max)
  sign_vec = sign(sapply(1:length(beta_max), function(x) coef1[x, beta_max_idx[x]]))

  scores = beta_max * sign_vec
  sources = colnames(mor_mat)

  tibble(score=scores, source=sources)
}
