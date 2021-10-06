#' Gene Set Variation Analysis (GSVA)
#'
#' @description
#' Calculates regulatory activities using GSVA.
#'
#' @details
#' This function is a wrapper for the method [GSVA::gsva()].
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param verbose Gives information about each calculation step. Default: FALSE.
#' @param method Method to employ in the estimation of gene-set enrichment.
#' scores per sample. By default this is set to gsva (Hänzelmann et al, 2013).
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
                     ...) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    regulons <- network %>%
        convert_to_gsva({{ .source }}, {{ .target }})

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
