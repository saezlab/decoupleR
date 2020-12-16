#' This function takes the elements of the `activity` column and calculates
#' precision-recall and ROC curves (depending on `curve`).
#' The `activity` column is populated with the output for each stat method and
#' results from the `run_benchmark()` function. Each of the elements
#' in `activity` are the results from runs of the `decouple()` wrapper.
#'
#' @param df run_benchmark roc column provided as input
#' @param downsampling logical flag indicating if the number of TN should be
#'   downsampled to the number of TP
#' @param times integer showing the number of downsampling
#' @param ranked logical flag indicating if input is derived from composite
#'   ranking that already took up-/downregulation (sign) into account
#' @param curve whether to return a Precision-Recall Curve ("PR") or ROC ("ROC")
#'
#' @return tidy data frame with precision, recall, auc, n, tp, tn and coverage
#'
#' @import yardstick
calc_curve = function(df, downsampling = T,
                         times = 1000,
                         ranked = F,
                         curve="PR") {

  if(curve=="PR"){
    res_col_1 <- "precision"
    res_col_2 <- "recall"
    curve_fun = yardstick::pr_curve
    auc_fun = yardstick::pr_auc
  }
  else if(curve=="ROC"){
    res_col_1 <- "sensitivity"
    res_col_2 <- "specificity"
    curve_fun = yardstick::roc_curve
    auc_fun = yardstick::roc_auc
  }


  if (ranked == T) {
    df = df %>% prepare_for_roc(., filter_tn = T, ranked = T)
  } else {
    df = df %>%
      prepare_for_roc(., filter_tn = T, ranked = F)
  }

  if (length(which(df$response == 0)) == nrow(df)){
    return(as_tibble(NULL))
  }

  tn = df %>% filter(response == 0)
  tp = df %>% filter(response == 1)

  feature_coverage = length(unique(df$tf))

  if (downsampling == T) {
    num_tp = nrow(tp)

    res = map_df(seq(from=1, to=times, by=1), function(i) {
      df_sub = sample_n(tn, num_tp, replace=TRUE) %>%
        bind_rows(tp)

      r_sub = df_sub %>%
        curve_fun(response, predictor)

      auc = df_sub %>%
        auc_fun(response, predictor) %>%
        pull(.estimate)

      res_sub = tibble({{ res_col_1 }} := r_sub %>% pull(res_col_1),
                       {{ res_col_2 }} := r_sub %>% pull(res_col_2),
                       th = r_sub$.threshold,
                       auc = auc,
                       n = length(which(df$response == 1)),
                       tp = nrow(tp),
                       tn = nrow(tn),
                       coverage = feature_coverage) %>%
        mutate_("run" = i)

    })
    # Get Average AUC
    res$auc <- sum(res$auc)/length(res$auc)
    res$tn <- nrow(tp)

  } else {
    r = df %>%
      curve_fun(response, predictor)
    auc = df %>%
      auc_fun(response, predictor)

    res = tibble({{ res_col_1 }} := r %>% pull(res_col_1),
                 {{ res_col_2 }} := r %>% pull(res_col_2),
                 th = r$.threshold,
                 auc = auc$.estimate,
                 n = length(which(df$response == 1)),
                 tp = nrow(tp),
                 tn = nrow(tn),
                 coverage = feature_coverage) %>%
      arrange(!!res_col_1, !!res_col_2)
  }

  return(res)
}


#' This function prepares each `activity` element for ROC curve calculation.
#'
#' @param df `activity` column elements - i.e. `decouple()` output.
#' @param filter_tn logical flag indicating if unnecessary true negatives should
#' be filtered out (unnecessary means that there are no true positives for a
#'   given source)
#' @param ranked logical flag indicating if input is derived from composite
#'   ranking that already took up-/downregulation (sign) into account
#'
#' @return tidy data frame with meta information for each experiment and the
#'   response and the predictor value which are required for roc curve analysis
#' @aliases source_set_columns
#'
prepare_for_roc = function(df, filter_tn = F, ranked = F) {
  res = df %>%
    dplyr::mutate(response = case_when(tf == target ~ 1,
                                       tf != target ~ 0),
                  predictor =  case_when(ranked == F ~ score*sign,
                                         ranked == T ~ score))
  res$response = factor(res$response, levels = c(1, 0))

  if (filter_tn == TRUE) {
    # Only TF which are perturbed and predicted are considered
    z = intersect(res$tf, res$target)
    res = res %>%
      filter(tf %in% z, target %in% z)
  }
  res %>%
    select(c(tf, id, response, predictor))
}
