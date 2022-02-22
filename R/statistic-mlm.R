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
#' MLM also provides the option to fit a regularized multivariate linear model
#' (ridge or lasso). Using the functions in the "glmnet" package, different
#' stringency of regularization (lambda) results in an ensemble of many solutions.
#' At each solution, the coefficients of regulators (sources) are scaled by the
#' absolute sum of all coefficients, representing the "importance" of each regulator.
#' For each regulator, its optimal "importance" (strongly positive or negative)
#' among the ensemble is used to represent its activity.

#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Deprecated parameter.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param minsize Integer indicating the minimum number of targets per source.
#' @param regularization Indicate whether apply regularization. NULL if no regularization. Set to 'ridge' or 'lasso' to use the corresponding regularization scheme implemented in glmnet. Can also set to numerical values between 0 (ridge) to 1 (lasso), and any value in between represent elastic net.
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
run_mlm <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    .mor = .data$mor,
                    .likelihood = .data$likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    minsize = 5,
                    regularization = NULL) {
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  alpha = NULL

  # Check the value in regularization is valid
  if (is.numeric(regularization)) {
    alpha = regularization
    stopifnot((alpha >= 0 ) && (alpha <= 1))
  }

  if (is.character(regularization)) {
    if (regularization == 'lasso') {
      alpha = 1
    } else if (regularization == 'ridge') {
      alpha = 0
    }
    else {
      stop("Invalid value for 'regularization")
    }
  }

  # Before to start ---------------------------------------------------------
  # Convert to standard tibble: source-target-mor.
  network <- network %>%
    rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
  network <- filt_minsize(rownames(mat), network, minsize)

  # Preprocessing -----------------------------------------------------------
  .fit_preprocessing(network, mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
  {
    .mlm_analysis(.$mat, .$mor_mat, alpha)
  } %>%
    ungroup()
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
.mlm_analysis <- function(mat, mor_mat, alpha) {
  mlm_evaluate_model <- partial(
    .f = .mlm_evaluate_model,
    mat = mat,
    mor_mat = mor_mat,
    alpha = alpha
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
.mlm_evaluate_model <- function(condition, mat, mor_mat, alpha=NULL) {
  if (is.null(alpha)) {
    fit <- lm(mat[ , condition] ~ mor_mat) %>%
      summary()
    scores <- as.vector(fit$coefficients[,3][-1])
    pvals <- as.vector(fit$coefficients[,4][-1])
    sources <- colnames(mor_mat)
    diff_n <- length(sources) - length(scores)
    if (diff_n > 0) {
      stop(stringr::str_glue('After intersecting mat and network, at least {diff_n} sources in the network are colinear with other sources.
      Cannot fit a linear model with colinear covariables, please remove them.
      Please run decoupleR::check_corr to see what regulators are correlated.
      Anything above 0.5 correlation should be removed.'))
    }
    tibble(score=scores, p_value=pvals, source=sources)
  }
  else {
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

    norm_beta = scale(abs(coef1)**(2-alpha), center=FALSE,
                      scale=colSums(abs(coef1)**(2-alpha)))
    beta_max = apply(norm_beta, 1, max)

    beta_max_idx = apply(norm_beta, 1, which.max)
    sign_vec = sign(sapply(1:length(beta_max), function(x) coef1[x, beta_max_idx[x]]))

    scores = beta_max * sign_vec
    sources = colnames(mor_mat)

    tibble(score=scores, p_value=NA, source=sources)
  }

}
