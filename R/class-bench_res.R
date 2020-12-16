#'  S4 Class used to format benchmark output.
#'
#' @field bench_res Formatted or non-formatted Benchmark output
#' @field summary Summary returned by the bench_sumplot functions - it contains
#' a summary table, roc plot, precision-recall curve (prc) plot, auroc heatmap,
#' and precision-recall area under the curve heatmap
#' @field design The input design tibble used to generate the benchmark results
#'
#' @exportClass BenchResult
setClass("BenchResult",
         slots=list(bench_res="tbl_df",
                    summary="list",
                    design="tbl_df"))
