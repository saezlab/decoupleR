#' Function to generate a consensus score between methods from the
#' result of decouple
#'
#' @param df `decouple` data frame result
#' @param include_time Should the time per statistic evaluated be informed?
#'
#' @return Updated tibble with the computed consensus score between methods
#'
#' @import purrr
#' @import RobustRankAggreg
#' @export
run_consensus <- function(df,
                          include_time=FALSE
                          ){
  start_time <- Sys.time()
  # Split df by samples
  cond_names <- unique(df$condition)
  lst_conds <- df %>%
    group_by(.data$condition) %>%
    group_split()
  names(lst_conds) <- cond_names

  # Split each sample by method
  run_id <- max(df$run_id)
  consensus <- lst_conds %>%
    # Generate a sorted list of sources per method
    map(function(df){
      df %>%
        group_by(.data$statistic) %>%
        group_split() %>%
        map(function(df){
          df %>%
            arrange(desc(abs(.data$score))) %>%
            select(.data$source) %>%
            pull()
        })
    }) %>%
    # Compute ranks
    map(function(lst){
      RobustRankAggreg::rankMatrix(lst) %>%
        RobustRankAggreg::aggregateRanks(rmat = .)
    }) %>%
    # Transform back to tibble
    map2(., names(.), function(df, cond){
      as_tibble(df) %>%
        rename('source' = .data$Name, 'p_value' = .data$Score) %>%
        mutate(score= -log10(.data$p_value),
               statistic = 'consensus',
               condition = cond,
               run_id = run_id + 1
        )
    }) %>%
    bind_rows()

  if (include_time) {
    consensus <- consensus %>%
      add_column(
        statistic_time = difftime(Sys.time(), start_time),
        .after = "score"
      )
  }

  # Join results
  result <- list(df, consensus) %>% bind_rows()

  result
}

