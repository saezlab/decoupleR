#' Function to format benchmarking results
#'
#' @param bench_res benchmarking results
#' @param .silent bool whether to silence warnings or not
#' @returns formatted benchmarking results
#' @importFrom rlang .data
#' @importFrom stringr str_glue_data
bench_format <- function(bench_res, .silent) {
  res_format <- bench_res %>%
    unnest(activity) %>%
    # convert filter_criteria from character to string
    rowwise() %>%
    mutate(filter_crit = paste0(unlist(filter_crit), collapse = "")) %>%
    ungroup() %>%
    # get statistic name
    mutate(statistic = activity %>%
             map(function(tib)
               unique(tib[["statistic"]]))) %>%
    unnest(statistic) %>%
    select(set_name, bench_name, filter_crit, statistic, activity)

  # Check and filter infinite values
  inf_sums <- lapply(res_format$activity,
                     function(x)
                       sum(is.infinite(x$score))) %>%
    setNames(
      paste(
        res_format$set_name,
        res_format$bench_name,
        res_format$statistic,
        sep = "_"
      )
    ) %>%
    enframe() %>% unnest(value)

  if (sum(inf_sums$value)) {
    res_format <- res_format %>%
      mutate(activity = activity %>%
               map(function(tib)
                 tib %>%
                   mutate_at(
                     vars(score), ~ replace(., is.infinite(.), 0)
                   )))

    if (!.silent)
      warning(
        inf_sums %>%
          filter(value > 0) %>%
          str_glue_data("{.$value} infinite values were filtered",
                        " in {.$name}. \n ")
      )

  }
  return(res_format)
}
