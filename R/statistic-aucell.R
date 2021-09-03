#' AUCell wrapper
#'
#' This function is a convenient wrapper for the workflow that calculates AUC
#' for each regulon in each cell. The workflow consists of the AUCell::AUCell_buildRankings
#' and AUCell::AUCell_calcAUC functions.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param aucMaxRank Threshold to calculate the AUC.
#' @param nCores Number of cores to use for computation.
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
#' run_aucell(mat, network)
run_aucell <- function(mat,
                       network,
                       .source = .data$tf,
                       .target = .data$target,
                       aucMaxRank = ceiling(0.05 * nrow(rankings)),
                       nCores = 1) {
  # Before to start ---------------------------------------------------------

  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  network <- network %>%
    convert_to_aucell({{ .source }}, {{ .target }})

  # Convert to absolute values
  mat <- abs(mat)

  # Analysis ----------------------------------------------------------------
  rankings <- exec(.fn = AUCell::AUCell_buildRankings,
                   exprMat = mat,
                   plotStats = FALSE,
                   verbose = FALSE,
                   nCores = nCores)

  exec(.fn = AUCell::AUCell_calcAUC,
       geneSets = network,
       rankings = rankings,
       verbose = FALSE,
       aucMaxRank = aucMaxRank,
       nCores = nCores
  ) %>%
    .extract_assay_auc() %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-tf ,names_to = "condition", values_to = "score") %>%
    add_column(statistic = "aucell", .before = 1)

}

.extract_assay_auc <- function(.a){
  SummarizedExperiment::assay(.a)
}
