#' Function that provides summary and plots for the benchmark run
#'
#' @param .res_tibble formatted bench result tibble with roc and prc columns
#' roc column: Reciever Operator Curve Results (calculated with yardstick)
#' prc column: Precision-Recall Curve Results (calculated with yardstick)
#' @return A summary list with TF coverage, ROC, AUROC, PRAUC, Run time,
#' ROC plots, and Heatmap plots
#' @import ggplot2 pheatmap
#' @importFrom rlang .data
get_bench_summary <- function(.res_tibble) {
  # get roc results
  roc <- format_roc(.res_tibble, "roc")

  # get PR roc results
  pr <- format_roc(.res_tibble, "prc")

  # Plot ROC
  roc_plot <- ggplot(roc, aes(x = 1-specificity,
                              y = sensitivity,
                              colour = run_key)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("FPR (1-specificity)") +
    ylab("TPR (sensitivity)")

  # Plot PR ROC
  pr_plot <- ggplot(pr, aes(x = recall, y = precision , colour = run_key)) +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    xlab("Recall/Sensitivity") +
    ylab("Precision")

  # Extract AUROC
  auroc_tibble <- .res_tibble %>%
    unnest(roc) %>%
    select(set_name, bench_name, filter_crit, statistic, auc) %>%
    distinct()

  # Plot AUROC
  auroc_plot <- auroc_tibble %>%
    unite("run_key", set_name, bench_name,
          statistic, filter_crit, remove = F) %>%
    ggplot(., aes(x = reorder(run_key, auc),
                  y = auc,
                  fill = run_key)) +
    geom_bar(stat = "identity") +
    xlab("networks") +
    ylab("AUROC") +
    coord_flip(ylim = c(0.5, 0.8)) +
    theme(legend.position = "none")

  # AUROC Heatmap
  auroc_heat <- auroc_tibble %>% get_auroc_heat()

  # Extract AU PRROC
  prauc_tibble <- .res_tibble %>%
    unnest(prc) %>%
    select(set_name, bench_name, filter_crit, statistic, auc) %>%
    distinct()

  # AU PR Heatmap
  pr_heat <- prauc_tibble %>% get_auroc_heat()

  # get computational time info
  comp_time <- .res_tibble %>%
    # get statistic time from activity
    mutate(statistic_time = activity %>%
             map(function(tib)
               tib %>%
                 select(statistic_time) %>%
                 unique)) %>%
    unnest(statistic_time) %>%
    # calculate regulon size
    group_by(set_name, bench_name, filter_crit) %>%
    mutate(regulon_time = sum(statistic_time)) %>%
    ungroup %>%
    select(set_name, bench_name, statistic,
           filter_crit, statistic_time, regulon_time)

  # Join AUROC, PRAUC, Coverage, and Comp time
  roc_cov <- roc %>%
                group_by(name_lvl) %>%
                summarise(source_cov = coverage,
                          condition_cov = n,
                          roc_neg=tn) %>%
                distinct() %>%
                ungroup() %>%
                separate(col="name_lvl",
                         into=c("set_name", "bench_name", "filter_crit"),
                         sep="\\.")

  pr_cov <- pr %>%
               group_by(name_lvl) %>%
               summarise(pr_neg=tn) %>%
               distinct() %>%
               ungroup() %>%
               separate(col="name_lvl",
                        into=c("set_name", "bench_name", "filter_crit"),
                        sep="\\.")


  summary_table <- auroc_tibble %>%
    inner_join(prauc_tibble %>%
                 rename(pr_auc = auc),
               by = c("set_name", "bench_name", "statistic", "filter_crit")) %>%
    inner_join(x=.,
               y=roc_cov,
               by = c("set_name", "bench_name", "filter_crit")) %>%
    inner_join(x=.,
               y=pr_cov,
               by = c("set_name", "bench_name", "filter_crit")) %>%
    inner_join(x=.,
               y=comp_time,
               by = c("set_name", "bench_name", "filter_crit", "statistic")) %>%
    distinct()

  bench_summary <- list(summary_table, roc_plot, pr_plot,
                        auroc_plot, auroc_heat, pr_heat)

  names(bench_summary) <- c("summary_table", "roc_plot", "pr_plot",
                            "auroc_plot", "auroc_heat", "pr_heat")

  return(bench_summary)
}



#' Helper function to format (Precision-Recall) Receiver Operator Curve results
#' @param .res_tibble formatted bench result tibble with added AUROC column
#' @param roc_column PR/ROC column to format
#' @return returns
format_roc <- function(.res_tibble, roc_column){
  .res_tibble %>%
    select(set_name, bench_name, filter_crit, statistic, roc_column) %>%
    unnest(roc_column) %>%
          unite("name_lvl", set_name, bench_name,
                filter_crit, remove = FALSE, sep = ".") %>%
          unite("run_key", set_name, bench_name,
                statistic, filter_crit, remove = FALSE)
}


#' Helper function to produce AUROC heatmap
#' @param auroc_tibble Tibble with calculated AUROC
#' @return returns an AUROC or Precision-Recall AUC heatmap
#' @import pheatmap ggplot2
get_auroc_heat <- function(auroc_tibble){
  auroc_tibble %>%
    select(statistic, auc, filter_crit, set_name, bench_name) %>%
    unite("name_lvl", set_name, bench_name, filter_crit) %>%
    pivot_wider(names_from = name_lvl, values_from = auc) %>%
    column_to_rownames(var = "statistic")  %>%
    pheatmap(.,
             cluster_rows = FALSE,
             treeheight_col = 0,
             treeheight_row = 0,
             display_numbers = TRUE,
             silent = TRUE,
             cluster_cols = FALSE)
}
