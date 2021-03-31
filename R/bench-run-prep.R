# Helper functions for run_benchmark

#' Helper Function to to generate the bools used to check if the current
#' locations/rds objects are the same as the previous one.
#' @param .design input tibble used to provide the experimental design for each
#' benchmark run
#' @keywords internal
#' @details This is used to limit the number of times that any of the
#' prerequsites is loaded.
format_design <- function(.design){
  .design %>%
    mutate(.source_bln = .data$source_loc %>% check_preced(),
           .expr_bln = .data$bexpr_loc %>% check_preced(),
           .meta_bln = .data$bmeta_loc %>% check_preced())
}


#' Helper Function that checks if the preceding vector element is the same
#' as the current element
#'
#' @param vector_loc character vector with directory paths
#'
#' @return logical values describing whether the location of the loaded files
#' has changes
#'
#' @keywords internal
check_preced <- function(vector_loc){
  tib_loc <- tibble(current=vector_loc, behind=lag(vector_loc))

  pmap_lgl(tib_loc, function(behind, current){
    ifelse(is.na(behind) || behind!=current, FALSE, TRUE)
  })
}


#' Helper Function to filter and format the gene set resource
#'
#' @param set_source Set Source (e.g. TF regulon sets, GO:term sets, etc)
#' @inheritParams input_tibble
#' @param .minsize minimum size of each set
#' @param .silent bool whether to silence wanring messages
#'
#' @importFrom stringr str_glue
#' @return returns a filtered and formatted set source
#' @details Filtering can be omitted if `filter_col` is `NA`.
filter_sets <- function(set_source,
                        source_col,
                        filter_col,
                        filter_crit,
                        .minsize,
                        .silent){
  n_duprows <- sum(duplicated(set_source))
  na_bool <- is.na(filter_col)
  print(na_bool)

  gs_filtered <- set_source %>%
    {
      if(na_bool){distinct(.)}
      else if(!na_bool){
        filter(., .data[[filter_col]] %in% filter_crit) %>%
          distinct_at(vars(-.data[[filter_col]]), .keep_all = FALSE)
        }
      } %>%
    group_by(.data[[source_col]]) %>%
    add_count() %>%
    filter(n >= .minsize) %>%
    ungroup()

  if (n_duprows & !.silent){
    warning(str_glue("{n_duprows} rows were duplicated in the set resource! ",
                     "{sum(duplicated(gs_filtered))} duplicated rows ",
                     "remain after filtering."))
  }
  return(gs_filtered)
}


#' `base::readRDS` helper function that enables loading files from urls
#' @inheritParams base::readRDS
#' @inheritDotParams base::readRDS
#' @param .url_bool bool whether the location is a url or not
#' @export
readRDS_helper <- function(file, .url_bool=FALSE, ...){
  if(.url_bool){
    readRDS(url(file, "rb", ...))
  } else{
    readRDS(file, ...)
  }
}
