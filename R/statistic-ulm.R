#' Univariate Linear Model (ULM)
#'
#' @description
#' Calculates regulatory activities using ULM.
#'
#' @details
#' ULM fits a linear model for each sample and regulator, where the observed
#' molecular readouts in mat are the response variable and the regulator weights
#' in net are the explanatory one. Target features with no associated weight
#' are set to zero. The obtained t-value from the fitted model is the activity
#' `ulm` of a given regulator.
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
#' run_ulm(mat, network, .source='tf')
run_ulm <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    .mor = .data$mor,
                    .likelihood = .data$likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
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
            .ulm_analysis(.$mat, .$mor_mat)
        } %>%
        ungroup()
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
        mutate(model = list(ulm_evaluate_model(.data$source, .data$condition)), 
               statistic='ulm', source=.data$source) %>%
        unnest(.data$model) %>%
        select(.data$statistic, .data$source, .data$condition, .data$score, .data$p_value)
}

#' Wrapper to run ulm one source (source) per sample (condition) at time
#'
#' @keywords internal
#' @noRd
.ulm_evaluate_model <- function(source, condition, mat, mor_mat) {
    #data <- cbind(data.frame(y=mat[ , condition]), mor_mat[, source])
    fit <- lm(mat[ , condition] ~ mor_mat[, source]) %>%
        summary()
    scores <- as.vector(fit$coefficients[,3][-1])
    pvals <- as.vector(fit$coefficients[,4][-1])
    tibble(score=scores, p_value=pvals)
}
