#' Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER)
#'
#' @description
#' Calculates regulatory activities using VIPER.
#'
#' @details
#' This function is a wrapper for the method `viper::viper()`.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @inheritDotParams viper::viper -eset -regulon -minsize
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import tibble
#' @import purrr
#' @import tidyr
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_viper(mat, network, .source='tf', verbose = FALSE)
run_viper <- function(mat,
                      network,
                      .source = .data$source,
                      .target = .data$target,
                      .mor = .data$mor,
                      .likelihood = .data$likelihood,
                      verbose = FALSE,
                      minsize = 0,
                      pleiotropy = T,
                      eset.filter = F,
                      ...) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    network <- network %>%
        convert_to_viper({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

    # Analysis ----------------------------------------------------------------
    exec(
        .fn = viper::viper,
        eset = mat,
        regulon = network,
        verbose = verbose,
        minsize = minsize,
        pleiotropy = pleiotropy,
        eset.filter = eset.filter,
        !!!list(...)
    ) %>%
        as.data.frame() %>%
        rownames_to_column("source") %>%
        pivot_longer(-.data$source, names_to = "condition", values_to = "score") %>%
        add_column(statistic = "viper", .before = 1) %>%
        mutate(p_value = 2*pnorm(-abs(score)))
}
