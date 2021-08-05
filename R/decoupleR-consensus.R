#' Function to generate a consensus score between methods from the
#' result of decouple
#'
#' @param df deocuple results
#' @param condition Column name for sample names
#' @param statistic Column name for statistic names
#'
#' @return Updated tibble with the computed consensus score bewteen methods
#'
#' @import purrr
#' @import RobustRankAggreg
#' @export
run_consensus <- function(df,
                          condition='condition',
                          statistic='statistic'){
  # Split df by samples
  cond_names <- unique(df$condition)
  lst_conds <- df %>%
    group_by(condition) %>%
    group_split()
  names(lst_conds) <- cond_names

  # Split each sample by method
  stats_names <- unique(df[[statistic]])
  consensus <- lst_conds %>%
    # Generate a sorted list of TFs per method
    map(function(df){
      df %>%
        group_by(statistic) %>%
        group_split() %>%
        map(function(df){
          df %>%
            arrange(desc(abs(score))) %>%
            select(tf) %>%
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
        rename('tf' = Name, 'score' = Score) %>%
        mutate(score= -log10(score),
               statistic = 'consensus',
               condition = cond,
               run_id = as.character(length(stats_names) + 1)
        )
    }) %>%
    bind_rows()

  # Join results
  result <- list(df, consensus) %>% bind_rows()
  result
}

