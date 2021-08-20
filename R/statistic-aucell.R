#' AUCell wrapper
#'
#' This function is a convenient wrapper for the AUCell workflow that consists of [AUCell::AUCell_buildRankings] and  [AUCell::AUCell_buildRankings] functions.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format


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
                       .target = .data$target) {
  # Check for NAs/Infs in mat
  check_nas_infs(mat)
  
  # Before to start ---------------------------------------------------------
  network <- network %>%
    convert_to_aucell({{ .source }}, {{ .target }})
  
  # Analysis ----------------------------------------------------------------
  rankings <- exec(.fn = AUCell::AUCell_buildRankings,
                   exprMat = mat,
                   plotStats = FALSE,
                   verbose = FALSE)
  
  
  exec(.fn = AUCell::AUCell_calcAUC,
       geneSets = network,
       rankings = rankings,
       verbose = FALSE
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