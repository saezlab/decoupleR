#' Benchmark input tibble containing the (experimental) design for each of the
#' benchmark runs corresponding to rows
#' @name input_tibble
#'
#' @details A tibble with locations, options, and filter options for
#' the desired benchmark setting
#'
#' @param set_name user-defined name of the set resource
#' @param bench_name user-defined name of the benchmark data
#' @param stats_list List of statistics to run
#' @param opts_list Named list containing the options for each stat. method
#' @param bexpr_loc benchmark expression data location (.rds format tibble)
#' @param bmeta_loc benchmark metadata location (.rds format tibble)
#' @param source_loc set source (e.g. network resource, gene ontology sets,
#'  kinase sets, etc.) location (.rds format tibble)
#' @param source_col name of the column with the source for the set source
#' @param target_col name of the column with the targets for the set source
#' @param filter_col name of the column by which we wish to filter
#' @param filter_crit criteria by which we wish to filter the filter column
NULL
