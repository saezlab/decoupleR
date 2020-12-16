#' Benchmark pipeline built on the statistical method wrapper decouple.
#'
#' @inheritParams format_design
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate roc and performance summary
#' @inheritParams filter_sets
#' @param .downsample_pr whether to downsample precision recall curve TNs
#' @param .downsample_roc whether to downsample ROC true negatives
#' @param .downsample_times downsampling iterattions
#' @inheritParams readRDS_helper
#'
#' @export
#' @importFrom rlang .data
#' @importFrom methods new
#' @importFrom stats reorder setNames
#' @import tibble tidyr dplyr tidyselect
#' @seealso See \link{input_tibble} for a description of the params/column
#' of .design (i.e. input tibble).
#' @return An S4 object of class BenchResult  \link{BenchResult}
run_benchmark <- function(.design,
                          .form = TRUE,
                          .perform = TRUE,
                          .minsize = 10,
                          .silent = TRUE,
                          .downsample_pr = FALSE,
                          .downsample_roc = FALSE,
                          .downsample_times = 100,
                          .url_bool = FALSE
                          ){
  res <- .design %>%
    format_design() %>%
    mutate(activity = pmap(.,
                           .f=function(set_name, bench_name,
                                       stats_list, opts_list,
                                       bexpr_loc, bmeta_loc, source_loc,
                                       source_col, target_col,
                                       filter_col, filter_crit,
                                       .source_bln, .expr_bln, .meta_bln){

      # Check_prereq
       if(!.expr_bln){
         .GlobalEnv$gene_expression <- readRDS_helper(bexpr_loc, .url_bool) %>%
           as.matrix()
       }
       if(!.meta_bln){
         .GlobalEnv$meta_data <- readRDS_helper(bmeta_loc, .url_bool)
       }
       if(!.source_bln){
         .GlobalEnv$set_source <- check_prereq(source_loc, source_col,
                                               filter_col, target_col,
                                               .url_bool)
       }

      # Filter set_source/network
      ss_filtered <- filter_sets(set_source, source_col,
                                 filter_col, filter_crit,
                                 .minsize, .silent)

      # Print Current Row/Run
      if(!.silent){
        .curr_row <- paste(set_name, bench_name,
                           paste0(unlist(filter_crit), collapse=""),
                           sep="_")
        message(str_glue("Currently Running: {.curr_row}"))
      }

      # Obtain Activity with decouple and format
      decouple(mat = gene_expression, network = ss_filtered,
               .source = source_col, .target = all_of(target_col),
               statistics = stats_list,
               args = opts_list)  %>%
        dplyr::rename(id=condition) %>%
        inner_join(meta_data, by="id")  %>%
        group_split(statistic, .keep=TRUE) %>%
        as.list()
      })) %>% {
      if(.form & !.perform) bench_format(., .silent)
      else if(.form & .perform) bench_format(., .silent) %>%
        mutate(roc = activity %>%
                 map(~calc_curve(df=.x,
                                    downsampling=.downsample_roc,
                                    times=.downsample_times,
                                    curve="ROC")),
               prc = activity %>%
                 map(~calc_curve(df=.x,
                                    downsampling=.downsample_pr,
                                    times=.downsample_times,
                                    curve="PR")))
      else .
    }

  if(.form & .perform){
    bench_result <-new("BenchResult",
                       bench_res=res,
                       summary=res %>% get_bench_summary(),
                       design=.design)
  }
  else{
    bench_result <-new("BenchResult",
                       bench_res=res,
                       summary=list(NULL),
                       design=.design)
  }
  return(bench_result)
}
