#' Function to format benchmarking results
#'
#' @param .bench_res benchmarking results
#' @param .silent bool whether to silence warnings or not
#' @returns formatted benchmarking results
#' @importFrom rlang .data
#' @importFrom stringr str_glue_data
#'
#' @export
#' @details If infinite values are present in the results, this function will
#'   notify the user.
bench_format <- function(.bench_res, .silent) {
  res_format <- .bench_res %>%
    unnest(.data$activity) %>%
    # convert filter_criteria from character to string
    rowwise() %>%
    mutate(filter_crit = paste0(unlist(.data$filter_crit), collapse = "")) %>%
    ungroup() %>%
    # get statistic name
    mutate(statistic = .data$activity %>%
             map(function(tib)
               unique(tib[["statistic"]]))) %>%
    unnest(.data$statistic) %>%
    select(.data$set_name, .data$bench_name, .data$filter_crit,
           .data$statistic, .data$activity)

  inf_sums <- res_format$activity %>%
    map(function(x) sum(is.infinite(x$score))) %>%
    setNames(
      paste(
        res_format$set_name,
        res_format$bench_name,
        res_format$statistic,
        sep = "_"
      )) %>%
    enframe() %>% unnest(.data$value)

  if (sum(inf_sums$value)) {
    res_format <- res_format %>%
      mutate(activity = .data$activity %>%
               map(function(tib) tib %>%
                   mutate_at(vars(.data$score),
                             ~ replace(.data, is.infinite(.data), 0))))

    if (!.silent)
      warning(
        inf_sums %>%
          filter(.data$value > 0) %>%
          str_glue_data("{.$value} infinite values were filtered",
                        " in {.$name}. \n ")
      )

  }
  return(res_format)
}
