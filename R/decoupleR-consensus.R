#' Consensus score between methods
#' @description 
#' Function to generate a consensus score between methods from the
#' result of the `decouple` function.
#' 
#' @param df `decouple` data frame result
#' @param include_time Should the time per statistic evaluated be informed?
#' @param seed Deprecated parameter.
#'
#' @return Updated tibble with the computed consensus score between methods
#'
#' @import purrr
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' results <- decouple(
#'    mat = mat,
#'    network = net,
#'    .source = "source",
#'    .target = "target",
#'    statistics = c("wmean", "ulm"),
#'    args = list(
#'             wmean = list(.mor = "mor", .likelihood = "likelihood"),
#'             ulm = list(.mor = "mor", .likelihood = "likelihood")
#'         ),
#'    consensus_score = FALSE,
#'    minsize = 0
#'    )
#' run_consensus(results)
run_consensus <- function(df,
                          include_time=FALSE,
                          seed = NULL
                          ){

    # NSE vs. R CMD check workaround
    condition <- score <- source <- statistic <- NULL

  start_time <- Sys.time()
  
  # Filter Infs
  is_inf <- !is.finite(df$score)
  if (any(is_inf)) {
    warning("Infs detected in score, will be set to NAs. This might effect the final
            consensus score since they will be ignored.")
    df <- df %>%
      dplyr::filter(!is_inf)
  }
  
  run_id <- max(df$run_id)
  consensus <- df %>%
    dplyr::group_by(statistic, condition) %>%
    dplyr::group_split() %>%
    purrr::map(function(df){
      pos <- df %>%
        dplyr::filter(score > 0) %>%
        rbind(., dplyr::mutate(., score=-score)) %>%
        dplyr::mutate(score=score / sd(score)) %>%
        dplyr::filter(score > 0)
      neg <- df %>%
        dplyr::filter(score <= 0) %>%
        rbind(., dplyr::mutate(., score=-score)) %>%
        dplyr::mutate(score=score / sd(score)) %>%
        dplyr::filter(score <= 0)
      rbind(pos,neg)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(condition, source) %>%
    dplyr::summarize(score=mean(score), .groups = 'drop') %>%
    dplyr::mutate(p_value = 2*stats::pnorm(-abs(score))) %>%
    tibble::add_column(
      statistic = 'consensus',
      .before = 'source'
    ) %>%
    tibble::add_column(
      run_id = run_id + 1,
      .before = 'statistic'
    )

  if (include_time) {
    consensus <- consensus %>%
      tibble::add_column(
        statistic_time = difftime(Sys.time(), start_time),
        .after = "score"
      )
  }
  consensus
}

