#' GSVA wrapper
#'
#' This function is a convenient wrapper for the [GSVA::gsva()] function.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @inheritDotParams GSVA::gsva -expr -gset.idx.list
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
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
#' run_gsva(mat, network, tf, target, verbose = FALSE)
run_gsva <- function(
    mat,
    network,
    .source = .data$tf,
    .target = .data$target,
    ...) {
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
        !!!list(...)
    ) %>%
        as.data.frame() %>%
        rownames_to_column(var = "tf") %>%
        pivot_longer(
            cols = -.data$tf,
            names_to = "condition",
            values_to = "score"
        ) %>%
        transmute(
            statistic = "gsva",
            .data$tf, .data$condition, .data$score
        )
}
