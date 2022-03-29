#' Univariate Decision Tree (UDT)
#'
#' @description
#' Calculates regulatory activities by using UDT.
#'
#' @details
#' UDT fits a single regression decision tree for each sample and regulator,
#' where the observed molecular readouts in mat are the response variable and
#' the regulator weights in net are the explanatory one. Target features with
#' no associated weight are set to zero. The obtained feature importance from
#' the fitted model is the activity `udt` of a given regulator.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Deprecated parameter.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param min_n An integer for the minimum number of data points in a node that
#' are required for the node to be split further.
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
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
                    min_n = 20,
                    seed = 42,
                    minsize = 5
) {
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
    withr::with_seed(seed, {
      .udt_analysis(.$mat, .$mor_mat, min_n)
    })
  }
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
.udt_analysis <- function(mat, mor_mat, min_n) {
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
  data <- tibble(x = mat[, condition, drop=FALSE] , y = mor_mat[, source])
  score <- rpart::rpart(y~x, data, minsplit=min_n) %>% pluck("variable.importance")

  if (is.null(score)) {
    score <- 0
    names(score) <- source
  }
  score
}
