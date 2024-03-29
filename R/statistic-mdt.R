#' Multivariate Decision Trees (MDT)
#'
#' @description
#' Calculates regulatory activities using MDT.
#'
#' @details
#'
#' MDT fits a multivariate regression random forest for each sample, where the
#' observed molecular readouts in mat are the response variable and the
#' regulator weights in net are the covariates. Target features with no
#' associated weight are set to zero. The obtained feature importances from the
#' fitted model are the activities `mdt` of the regulators in net.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Deprecated parameter.
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
#' @importFrom magrittr %<>% %>%
#' @importFrom withr with_seed
#' @importFrom parallelly availableCores
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_mdt(mat, net, minsize=0)
run_mdt <- function(mat,
                    network,
                    .source = source,
                    .target = target,
                    .mor = mor,
                    .likelihood = likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    trees = 10,
                    min_n = 20,
                    nproc = availableCores(),
                    seed = 42,
                    minsize = 5
) {

    # NSE vs. R CMD check workaround
    condition <- likelihood <- mor <- score <- source <- target <- NULL

    # Check for NAs/Infs in mat
    mat %<>% check_nas_infs

    # Before to start ---------------------------------------------------------
    # Convert to standard tibble: source-target-mor.
    network %>%
    rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }}) %>%
    filt_minsize(rownames(mat), ., minsize) %>%
    # Preprocessing -----------------------------------------------------------
    .fit_preprocessing(mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
    {with_seed(seed, {.mdt_analysis(.$mat, .$mor_mat, trees, min_n, nproc)})}

}


#' Wrapper to execute run_mdt() logic one finished preprocessing of data
#'
#'
#' @inheritParams run_mdt
#' @param mor_mat
#'
#' @inherit run_mdt return
#' @importFrom purrr partial
#' @importFrom tidyr expand_grid
#' @importFrom dplyr rowwise reframe mutate arrange relocate
#' @keywords intern
#' @noRd
.mdt_analysis <- function(mat, mor_mat, trees, min_n, nproc) {

     # NSE vs R CMD check workaround
    condition <- NULL
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
    rowwise(condition) %>%
    reframe(
        score = mdt_evaluate_model(condition),
        source = colnames(mor_mat),
    ) %>%
    mutate(
        statistic = "mdt",
        source, condition, score,
        .before = 1L
    ) %>%
    arrange(source) %>%
    relocate(source, .after = 1L)

}


#' Wrapper to run mdt per a sample (condition) at time
#'
#' @importFrom purrr pluck
#' @importFrom magrittr %>%
#' @keywords internal
#' @noRd
.mdt_evaluate_model <- function(condition, mat, mor_mat, trees, min_n, nproc) {

    ranger::ranger(
        condition ~ .,
        data = data.frame(condition = mat[, condition], mor_mat),
        num.trees = trees,
        importance = "impurity",
        min.node.size = min_n,
        num.threads = nproc
    ) %>%
    pluck("variable.importance")

}
