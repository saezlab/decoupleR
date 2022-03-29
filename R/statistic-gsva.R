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
#' @param minsize Integer indicating the minimum number of targets per source.
#' @inheritDotParams GSVA::gsva -expr -gset.idx.list
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_gsva(mat, network, .source='tf', verbose = FALSE)
run_gsva <- function(mat,
                     network,
                     .source = .data$source,
                     .target = .data$target,
                     verbose = FALSE,
                     method = "gsva",
                     minsize = 5,
                     ...) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    network <- network %>%
        rename_net({{ .source }}, {{ .target }})
    network <- filt_minsize(rownames(mat), network, minsize)
    regulons <- extract_sets(network)

    # Analysis ----------------------------------------------------------------
    exec(
        .fn = GSVA::gsva,
        expr = mat,
        gset.idx.list = regulons,
        min.sz = 1,
        max.sz = Inf,
        verbose = verbose,
        method = method,
        !!!list(...)
    ) %>%
        as.data.frame() %>%
        rownames_to_column(var = "source") %>%
        pivot_longer(
            cols = -.data$source,
            names_to = "condition",
            values_to = "score"
        ) %>%
        transmute(
            statistic = "gsva",
            .data$source, .data$condition, .data$score
        )
}
