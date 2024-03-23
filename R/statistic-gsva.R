#' Gene Set Variation Analysis (GSVA)
#'
#' @description
#' Calculates regulatory activities using GSVA.
#'
#' @details
#' GSVA (Hänzelmann et al., 2013) starts by transforming the input molecular
#' readouts in mat to a readout-level statistic using Gaussian kernel estimation
#'  of the cumulative density function. Then, readout-level statistics are
#'  ranked per sample and normalized to up-weight the two tails of the rank
#'  distribution. Afterwards, an enrichment score `gsva` is calculated
#'  using a running sum statistic that is normalized by subtracting the largest
#'  negative estimate from the largest positive one.
#'
#'  Hänzelmann S. et al. (2013) GSVA: gene set variation analysis for microarray
#'   and RNA-seq data. BMC Bioinformatics, 14, 7.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param verbose Gives information about each calculation step. Default: FALSE.
#' @param method Method to employ in the estimation of gene-set enrichment.
#' scores per sample. By default this is set to gsva (Hänzelmann et al, 2013).
#' Further available methods are "plage", "ssgsea" and "zscore". Read more in
#' the manual of \code{\link{GSVA::gsva}}.
#' @param minsize Integer indicating the minimum number of targets per source.
#' Must be greater than 0.
#' @param maxsize Integer indicating the maximum number of targets per source.
#' @inheritDotParams GSVA::gsvaParam -exprData -geneSets -minSize -maxSize
#' @inheritDotParams GSVA::ssgseaParam -exprData -geneSets -minSize -maxSize
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @importFrom rlang !!! exec
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom tibble rownames_to_column
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_gsva(mat, net, minsize=1, verbose = FALSE)
run_gsva <- function(mat,
                     network,
                     .source = source,
                     .target = target,
                     verbose = FALSE,
                     method = c("gsva", "plage", "ssgsea", "zscore"),
                     minsize = 5L,
                     maxsize = Inf,
                     ...) {

    # NSE vs. R CMD check workaround
    condition <- score <- source <- target <- NULL

    if (minsize < 1L) {
        paste(
            'decoupleR::run_gsva: `minsize` must be greater than 0.',
            'Using 1 as minimum number of targets per source.'
        ) %>%
        warning(call. = FALSE)
        minsize <- 1L
    }

    param <- tryCatch(
        get(
            sprintf('%sParam', method[1L]),
            envir = asNamespace('GSVA'),
            inherits = FALSE
        ),
        error = function(e) {
            stop(sprintf(
                'No such method in GSVA: `%s`. To learn more check ?gsva.',
                method
            ))
        }
    )

    # Check for NAs/Infs in mat
    mat <- check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    network <- network %>%
        rename_net({{ .source }}, {{ .target }})
    network <- filt_minsize(rownames(mat), network, minsize)
    regulons <- extract_sets(network)

    # Analysis ----------------------------------------------------------------
    GSVA::gsva(
        expr = exec(
            param,
            exprData = mat,
            geneSets = regulons,
            minSize = minsize,
            maxSize = maxsize,
            !!!list(...)
        ),
        verbose = verbose
    ) %>%
    as.data.frame() %>%
    rownames_to_column(var = "source") %>%
    pivot_longer(
        cols = -source,
        names_to = "condition",
        values_to = "score"
    ) %>%
    mutate(
        statistic = "gsva",
        source, condition, score,
        .keep = "none",
        .before = 1L
    )

}
