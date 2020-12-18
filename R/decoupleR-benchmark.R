#' Benchmark pipeline built on the statistical method wrapper \link{decouple}.
#'
#' @inheritParams format_design
#' @param .form bool whether to format or not
#' @param .perform bool whether to calculate ROC and performance summary
#' @inheritParams filter_sets
#' @param .downsample_pr whether to downsample precision recall curve TNs
#' @param .downsample_roc whether to downsample ROC true negatives
#' @param .downsample_times downsampling iterations
#' @inheritParams readRDS_helper
#' @import tibble tidyr dplyr tidyselect
#' @seealso See \link{input_tibble} for a description of the params/columns
#'   of .design (i.e. input tibble).
#' @export
#' @importFrom rlang .data
#' @importFrom stats reorder setNames
#' @importFrom methods new
#' @return An S4 object of \link{BenchResult-class}
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

  bench_env <- new.env()

  res <- .design %>%
    format_design() %>%
    mutate(activity = pmap(.l=.,
                           .f=function(set_name, bench_name,
                                       stats_list, opts_list,
                                       bexpr_loc, bmeta_loc, source_loc,
                                       source_col, target_col,
                                       filter_col, filter_crit,
                                       .source_bln, .expr_bln, .meta_bln){

      # Check_prereq
       if(!.expr_bln){
         bench_env$gene_expression <- readRDS_helper(bexpr_loc, .url_bool) %>%
           as.matrix()
       }
       if(!.meta_bln){
         bench_env$meta_data <- readRDS_helper(bmeta_loc, .url_bool)
       }
       if(!.source_bln){
         bench_env$set_source <- readRDS_helper(source_loc, .url_bool)
       }

      # Filter set_source/network
      ss_filtered <- filter_sets(bench_env$set_source, source_col,
                                 filter_col, filter_crit,
                                 .minsize, .silent)

      # Show Current Row/Run
      if(!.silent){
        .curr_row <- paste(set_name, bench_name,
                           paste0(unlist(filter_crit), collapse=""),
                           sep="_")
        message(str_glue("Currently Running: {.curr_row}"))
      }

      # Obtain Activity with decouple and format
      decouple(mat = bench_env$gene_expression, network = ss_filtered,
               .source = source_col, .target = all_of(target_col),
               statistics = stats_list, args = opts_list)  %>%
        dplyr::rename(id=.data$condition) %>%
        inner_join(bench_env$meta_data, by="id")  %>%
        group_split(.data$statistic, .keep=TRUE) %>%
        as.list()
      })) %>% {
      if(.form & !.perform) bench_format(., .silent)
      else if(.form & .perform) bench_format(., .silent) %>%
        mutate(roc = .data$activity %>%
                 map(~calc_curve(df=.x,
                                    downsampling=.downsample_roc,
                                    times=.downsample_times,
                                    curve="ROC")),
               prc = .data$activity %>%
                 map(~calc_curve(df=.x,
                                    downsampling=.downsample_pr,
                                    times=.downsample_times,
                                    curve="PR")))
      else .
    }

  if(.form & .perform){
    bench_result <- methods::new("BenchResult",
                                 bench_res=res,
                                 summary=res %>% get_bench_summary(),
                                 design=.design)
  }
  else{
    bench_result <- methods::new("BenchResult",
                                 bench_res=res,
                                 summary=list(NULL),
                                 design=.design)
  }
  return(bench_result)
}
