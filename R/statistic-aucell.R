#' AUCell
#'
#' @description
#' Calculates regulatory activities using Area Under the Curve (AUC) from AUCell
#'
#' @details
#' This function is a wrapper for the method `AUCell`. It uses the
#' "Area Under the Curve" (AUC) to calculate whether a critical subset of input
#' molecular features is enriched for each sample.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param aucMaxRank Threshold to calculate the AUC.
#' @param nproc Number of cores to use for computation.
#'
#' @family decoupleR statistics
#' @export
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_aucell(mat, network, .source='tf')
run_aucell <- function(mat,
                       network,
                       .source = .data$source,
                       .target = .data$target,
                       aucMaxRank = ceiling(0.05 * nrow(rankings)),
                       nproc = 4,
                       seed = 42
) {
  # Before to start ---------------------------------------------------------
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  network <- network %>%
    convert_to_aucell({{ .source }}, {{ .target }})

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
