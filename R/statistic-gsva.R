#' GSVA wrapper
#'
#' This function is a convenient wrapper for the [GSVA::gsva()] function.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @inheritDotParams GSVA::gsva -expr -gset.idx.list
#'
#' @return
#'  A long format tibble of the enrichment scores for each tf across the conditions.
#'  Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#'  5. `statistic_time`: Internal execution time indicator.
#' @family decoupleR statistics
#' @export
run_gsva <- function(
    mat,
    network,
    .source = .data$tf,
    .target = .data$target,
    ...) {
    # Before to start ---------------------------------------------------------
    .start_time <- Sys.time()

    regulons <- network %>%
        convert_to_gsva({{ .source }}, {{ .target }})

    # Analysis ----------------------------------------------------------------
    exec(
        .fn = GSVA::gsva,
        expr = mat,
        gset.idx.list = regulons,
        !!!list(...)
    ) %>%
        as.data.frame() %>%
        rownames_to_column(var = "tf") %>%
        pivot_longer(cols = -.data$tf, names_to = "condition", values_to = "score") %>%
        transmute(statistic = "gsva", .data$tf, .data$condition, .data$score) %>%
        mutate(statistic_time = difftime(Sys.time(), .start_time))
}
