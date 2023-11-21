
    # NSE vs. R CMD check workaround
    c_score <- condition <- corr_wmean <- likelihood <- mor <- norm_wmean <-    null_distribution <- null_mean <- null_sd <- p_value <- score <- source <-    statistic <- target <- value <- weight <- wmean <- z_score <- NULL

#' AUCell
#'
#' @description
#' Calculates regulatory activities using AUCell.
#'
#' @details
#' AUCell (Aibar et al., 2017) uses the Area Under the Curve (AUC) to calculate
#' whether a set of targets is enriched within the molecular readouts of each
#' sample. To do so, AUCell first ranks the molecular features of each sample
#' from highest to lowest value, resolving ties randomly. Then, an AUC can be
#' calculated using by default the top 5% molecular features in the ranking.
#' Therefore, this metric, `aucell`, represents the proportion of
#' abundant molecular features in the target set, and their relative abundance
#' value compared to the other features within the sample.
#'
#' Aibar S. et al. (2017) Scenic: single-cell regulatory network inference and
#' clustering. Nat. Methods, 14, 1083–1086.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param aucMaxRank Threshold to calculate the AUC.
#' @param nproc Number of cores to use for computation.
#' @param seed A single value, interpreted as an integer, or NULL for random
#' number generation.
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom parallelly availableCores
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_aucell(mat, net, minsize=0, nproc=1, aucMaxRank=3)
run_aucell <- function(mat,
                       network,
                       .source = source,
                       .target = target,
                       aucMaxRank = ceiling(0.05 * nrow(rankings)),
                       nproc = availableCores(),
                       seed = 42,
                       minsize = 5
) {

    # NSE vs. R CMD check workaround
    source <- target <- NULL

  # Before to start ---------------------------------------------------------
  # Check for NAs/Infs in mat
  mat <- check_nas_infs(mat)

  network <- network %>%
    rename_net({{ .source }}, {{ .target }})
  network <- filt_minsize(rownames(mat), network, minsize)
  network <- extract_sets(network)

  # Convert to absolute values
  mat <- abs(mat)

  # Analysis ----------------------------------------------------------------
  withr::with_seed(seed, {
    rankings <- exec(.fn = AUCell::AUCell_buildRankings,
                     exprMat = mat,
                     plotStats = FALSE,
                     verbose = FALSE,
                     nCores = nproc)
  })

  withr::with_seed(seed, {
    exec(.fn = AUCell::AUCell_calcAUC,
         geneSets = network,
         rankings = rankings,
         verbose = FALSE,
         aucMaxRank = aucMaxRank,
         nCores = nproc
    )
  }) %>%
    .extract_assay_auc() %>%
    as.data.frame() %>%
    rownames_to_column("source") %>%
    pivot_longer(-source ,names_to = "condition", values_to = "score") %>%
    add_column(statistic = "aucell", .before = 1)

}

.extract_assay_auc <- function(.a){
  SummarizedExperiment::assay(.a)
}
