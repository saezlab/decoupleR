
    # NSE vs. R CMD check workaround
    ES <- NES <- condition <- p_value <- pathway <- pval <- score <- source <-    statistic <- target <- NULL

#' Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER)
#'
#' @description
#' Calculates regulatory activities using VIPER.
#'
#' @details
#' VIPER (Alvarez et al., 2016) estimates biological activities by performing a
#' three-tailed enrichment score calculation. For further information check the
#' supplementary information of the decoupler manuscript or the original
#' publication.
#' 
#' Alvarez M.J.et al. (2016) Functional characterization of somatic mutations
#' in cancer using network-based inference of protein activity. Nat. Genet.,
#' 48, 838–847.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param verbose Logical, whether progression messages should be printed in
#' the terminal.
#' @param pleiotropy Logical, whether correction for pleiotropic regulation
#' should be performed.
#' @param eset.filter Logical, whether the dataset should be limited only to
#' the genes represented in the interactome.
#' @param minsize Integer indicating the minimum number of targets per source.
#' @inheritDotParams viper::viper -eset -regulon -verbose -minsize -pleiotropy -eset.filter
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
#' @import tibble
#' @import purrr
#' @import tidyr
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_viper(mat, net, minsize=0, verbose = FALSE)
run_viper <- function(mat,
                      network,
                      .source = source,
                      .target = target,
                      .mor = mor,
                      .likelihood = likelihood,
                      verbose = FALSE,
                      minsize = 5,
                      pleiotropy = TRUE,
                      eset.filter = FALSE,
                      ...) {

    # NSE vs. R CMD check workaround
    likelihood <- mor <- score <- source <- target <- NULL

    # Check for NAs/Infs in mat
    mat <- check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    network <- network %>%
        rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
    network <- filt_minsize(rownames(mat), network, minsize)
    # Normalize mor between -1 and 1
    network <- network %>%
        dplyr::group_by(source) %>%
        dplyr::group_modify(function(.x, .y){
            n_max <- max(abs(.x$mor))
            .x$mor <- .x$mor / n_max
            .x
        })
    # Transform to viper format
    network <- network %>%
        dplyr::mutate(mor = mor) %>%
        split(.$source) %>%
        purrr::map(~ {
            list(
                tfmode = purrr::set_names(.x$mor, .x$target),
                likelihood = .x$likelihood
            )
        })

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
        pivot_longer(-source, names_to = "condition", values_to = "score") %>%
        add_column(statistic = "viper", .before = 1) %>%
        mutate(p_value = 2*stats::pnorm(-abs(score)))
}
