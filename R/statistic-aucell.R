#' AUCell wrapper
#'
#' This function is a convenient wrapper for the workflow that calculates AUC
#' for each regulon in each cell. The workflow consists of the AUCell::AUCell_buildRankings
#' and AUCell::AUCell_calcAUC functions.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param nCores Number of cores to use for computation.
#' @param weights Whether is network weighted
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
                       weights = FALSE,
                       nCores = 1) {
  # Check for NAs/Infs in mat
  check_nas_infs(mat)
  
  if (-1 %in% network$mor | weights){
    # Analysis ----------------------------------------------------------------
    network %>%
      convert_f_defaults(
        tf = {{ .source }},
        target = {{ .target }},
        mor = {{ .mor }},
        likelihood = {{ .likelihood }}
      ) %>%
      filter(target %in% rownames(mat)) %>% # Overlap between genes in the network and genes in the expression matrix
      group_by(tf) %>%
      group_modify(~ .one_TF_aucell_mor(.x, mat, nCores=nCores), .keep = TRUE) %>%
      pivot_longer(-tf , names_to = "condition",  values_to = "score") %>%
      add_column(statistic = "aucell", .before = 1) %>%
      ungroup()
    
  }
  else {
    # Before to start ---------------------------------------------------------
    network <- network %>%
      convert_to_aucell({{ .source }}, {{ .target }})
    
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
         nCores = nCores
    ) %>%
      .extract_assay_auc() %>%
      as.data.frame() %>%
      rownames_to_column("tf") %>%
      pivot_longer(-tf ,names_to = "condition", values_to = "score") %>%
      add_column(statistic = "aucell", .before = 1)
    
  }
}

.one_TF_aucell_mor <- function(network, mat, nCores){
  # Multiply mat by (likelihood*mor)
  mat[rownames(mat) %in% network$target,] <- (network$likelihood * network$mor) * mat[rownames(mat) %in% network$target,]
  
  # Calculate rankings
  .rankings <- AUCell::AUCell_buildRankings(mat, nCores=nCores, plotStats=FALSE, verbose=FALSE)
  
  # Convert network into named list for the AUCell_calcAUC function
  .network_aucell <- network %>%
    group_by(tf) %>%
    summarise(regulons = rlang::set_names(list(target), tf[1]),
              .groups = "drop"
    ) %>%
    pull(regulons)
  
  # Calculate AUC
  AUCell::AUCell_calcAUC(.network_aucell, .rankings, nCores=nCores, verbose=FALSE) %>%
    .extract_assay_auc() %>%
    as.data.frame()
  
}

.extract_assay_auc <- function(.a){
  SummarizedExperiment::assay(.a)
}