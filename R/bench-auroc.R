#' This function takes the elements of the `activity` column and calculates
#'   precision-recall and ROC curves (depending on `curve`).
#' The `activity` column is populated with the output for each stat method and
#'   results from the `run_benchmark()` function. Each of the elements
#'   in `activity` are results from runs of the \link{decouple} wrapper.
#'
#' @param df run_benchmark roc column provided as input
#' @param downsampling logical flag indicating if the number of Negatives
#'   should be downsampled to the number of Positives
#' @param times integer showing the number of downsampling
#' @param ranked logical flag indicating if input is derived from composite
#'   ranking that already took up-/downregulation (sign) into account
#' @param curve whether to return a Precision-Recall Curve ("PR") or ROC ("ROC")
#' @param seed An integer to set the RNG state for random number generation. Use
#'   NULL for random number generation.
#'
#' @return tidy data frame with precision, recall, auc, n, cp, cn and coverage
#'   in the case of PR curve; or sensitivity and specificity, auc, n, cp, cn
#'   and coverage in the case of ROC.
#' @import yardstick
calc_curve = function(df,
                      downsampling = FALSE,
                      times = 1000,
                      ranked = FALSE,
                      curve = "ROC",
                      seed = 420){

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

  if (sum(which(df$response == 0)) == nrow(df)){
    return(as_tibble(NULL))
  }

  cn = df %>% filter(response == 0)
  cp = df %>% filter(response == 1)

  feature_coverage = length(unique(df$tf))

  if (downsampling == TRUE) {
    num_tp = nrow(cp)

    res = map_df(seq(from=1, to=times, by=1), function(i) {
      set.seed(seed)
      df_sub = sample_n(cn, num_tp, replace=TRUE) %>%
        bind_rows(cp)

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
                       cp = nrow(cp),
                       cn = nrow(cn),
                       coverage = feature_coverage) %>%
        mutate("run" = i)

    })
    # Get Average AUC
    res$auc <- sum(res$auc)/length(res$auc)
    res$cn <- nrow(cp)

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
                 cp = nrow(cp),
                 cn = nrow(cn),
                 coverage = feature_coverage) %>%
      arrange(!!res_col_1, !!res_col_2)
  }

  return(res)
}


#' Helper function used to prepare `activity` elements or \link{decouple}
#' outputs for \link{calc_curve}. This is done by keeping only the the perturbed
#' or predicted `sources` (or TFs) and assigning the `score` (or statistical
#' method results e.g. Normalized Enrichment Score) as the `predictor`.
#'
#' @param df `activity` column elements - i.e. `decouple()` output.
#' @param filter_tn logical flag indicating if unnecessary negatives should
#' be filtered out
#' @param ranked logical flag indicating if input is derived from composite
#'   ranking that already took up-/downregulation (sign) into account
#'
#' @return tidy data frame with meta information for each experiment and the
#'   response and the predictor value which are required for ROC and
#'   PR curve analysis
prepare_for_roc = function(df, filter_tn = FALSE, ranked = FALSE) {
  res = df %>%
    dplyr::mutate(response = case_when(tf == target ~ 1,
                                       tf != target ~ 0),
                  predictor =  case_when(ranked == FALSE ~ score*sign,
                                         ranked == TRUE ~ score))
  res$response = factor(res$response, levels = c(1, 0))

  if (filter_tn == TRUE) {
    z = intersect(res$tf, res$target)
    res = res %>%
      filter(tf %in% z, target %in% z)
  }
  res %>%
    select(tf, id, response, predictor)
}
