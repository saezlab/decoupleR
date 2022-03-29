#' Multivariate Linear Model (MLM)
#'
#' @description
#' Calculates regulatory activities using MLM.
#'
#' @details
#' 
#' MLM fits a multivariate linear model for each sample, where the observed
#' molecular readouts in mat are the response variable and the regulator weights
#' in net are the covariates. Target features with no associated weight are set
#' to zero. The obtained t-values from the fitted model are the activities 
#' (`mlm`) of the regulators in net.
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
run_mlm <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    .mor = .data$mor,
                    .likelihood = .data$likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    minsize = 5) {
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
    .mlm_analysis(.$mat, .$mor_mat)
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
  fit <- lm(mat[ , condition] ~ mor_mat) %>%
    summary()
  scores <- as.vector(fit$coefficients[,3][-1])
  pvals <- as.vector(fit$coefficients[,4][-1])
  sources <- colnames(mor_mat)
  diff_n <- length(sources) - length(scores)
  if (diff_n > 0) {
    stop(stringr::str_glue('After intersecting mat and network, at least {diff_n} sources in the network are colinear with other sources.
      Cannot fit a linear model with colinear covariables, please remove them.
      Please run decoupleR::check_corr to see what regulators are correlated.'))
  }
  tibble(score=scores, p_value=pvals, source=sources)
}
