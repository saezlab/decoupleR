#' Fast Gene Set Enrichment Analysis (FGSEA)
#'
#' @description
#' Calculates regulatory activities using FGSEA.
#'
#' @details
#' This function is a wrapper for the method `fgsea::fgsea`.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param times How many permutations to do?
#' @param nproc Number of cores to use for computation.
#' @param seed A single value, interpreted as an integer, or NULL.
#' @inheritDotParams fgsea::fgseaMultilevel -pathways -stats -nPermSimple -nproc
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
#' run_fgsea(mat, network, .source='tf', nproc=1)
run_fgsea <- function(mat,
                      network,
                      .source = .data$source,
                      .target = .data$target,
                      times = 100,
                      nproc = 4,
                      seed = 42,
                      ...) {
  # Check for NAs/Infs in mat
  check_nas_infs(mat)

  regulons <- network %>%
    convert_to_fgsea({{ .source }}, {{ .target }})

  conditions <- colnames(mat) %>%
    set_names()

  map_dfr(.x = conditions, .f = ~ {
    stats <- mat[, .x]
    options <- list(
      pathways = regulons,
      stats = stats,
      nPermSimple = times,
      nproc = nproc
      )
    withr::with_seed(seed, {
      result <- suppressWarnings(do.call(what = fgsea::fgsea, args = options))
    })
  }, .id = "condition") %>%
    select(.data$pathway, .data$condition, .data$ES, .data$NES, .data$pval) %>%
    tidyr::pivot_longer(cols=c("ES","NES"), names_to ="statistic", values_to="score") %>%
    mutate(statistic=if_else(.data$statistic=='ES', 'fgsea', 'norm_fgsea')) %>%
    rename('source'=.data$pathway, 'p_value'=.data$pval) %>%
    select(.data$statistic, .data$source, .data$condition, .data$score, .data$p_value) %>%
    mutate(score = replace_na(.data$score, 0))
}
