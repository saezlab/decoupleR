#' Multivariate Linear Model (MLM)
#'
#' @description
#' Calculates regulatory activities using MLM.
#'
#' @details
#' 
#' MLM fits a multivariate linear model for each sample, where the observed
#' molecular readouts in mat are the response variable and the regulator weights
#' in net are the covariates. Target features with no associated weight are set
#' to zero. The obtained t-values from the fitted model are the activities 
#' (`mlm`) of the regulators in net.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Deprecated parameter.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#'
#' @import dplyr
#' @import purrr
#' @import tibble
#' @import tidyr
#' @importFrom stats coef lm summary.lm
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_mlm(mat, net, minsize=0)
run_mlm <- function(mat,
                    network,
                    .source = source,
                    .target = target,
                    .mor = mor,
                    .likelihood = likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    minsize = 5) {

    # NSE vs. R CMD check workaround
    condition <- likelihood <- mor <- p_value <- score <- source <- statistic    <- target <- NULL

  # Check for NAs/Infs in mat
  mat <- check_nas_infs(mat)
  
  # Before to start ---------------------------------------------------------
  # Convert to standard tibble: source-target-mor.
  network <- network %>%
    rename_net({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})
  network <- filt_minsize(rownames(mat), network, minsize)
  
  # Preprocessing -----------------------------------------------------------
  .fit_preprocessing(network, mat, center, na.rm, sparse) %>%
    # Model evaluation --------------------------------------------------------
  {
    .mlm_analysis(.$mat, .$mor_mat)
  } %>%
    ungroup()
}


#' Wrapper to execute run_mlm() logic one finished preprocessing of data
#'
#' Fit a linear regression between the value of expression and the profile of its targets.
#'
#' @inheritParams run_mlm
#' @param mor_mat
#'
#' @inherit run_mlm return
#' @keywords intern
#' @noRd


.mlm_analysis <- function(mat, mor_mat) {
  
    # run all linear models at the same time:
    res_all <- lm(mat ~ mor_mat) %>% summary()
    
    if(ncol(mat) == 1){
        # in case of a single condition, the summary of lm returns the table instead of 
        # list of tables.
        #
        res_all = list(res_all)
    }
    names(res_all) <- colnames(mat)
    
    # summary is a list for each condition. Get the info we need: 
    res_new <- res_all %>% lapply(X = ., function(fit){
        
        scores <- as.vector(fit$coefficients[,3][-1])
        pvals <- as.vector(fit$coefficients[,4][-1])
        sources <- colnames(mor_mat)
        diff_n <- length(sources) - length(scores)
        if (diff_n > 0) {
            stop(stringr::str_glue('After intersecting mat and network, at least {diff_n} sources in the network are colinear with other sources.
      Cannot fit a linear model with colinear covariables, please remove them.
      Please run decoupleR::check_corr to see what regulators are correlated.'))
        }
        tibble(score=scores, p_value=pvals, source=sources)
    }) %>% bind_rows(.id = "condition") %>%
        mutate(statistic = "mlm", .before= 1) %>%
        select(statistic, source, condition,
               score, p_value)
    return(res_new)
}




