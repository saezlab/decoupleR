#' z-score
#'
#' @description
#' Calculates regulatory activities using a z-score as descibed in KSEA or RoKAI.
#'
#' @details
#' The z-score calculates the mean of the molecular features of the known targets
#' for each regulator and adjusts it for the number of identified targets for the 
#' regulator, the standard deviation of all molecular features (RoKAI), as well as
#' the mean of all moleculare features (KSEA). 
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Deprecated parameter.
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param minsize Integer indicating the minimum number of targets per source.
#' @param flavor Whether the calculation should be based on RoKAI (default) or
#'  KSEA.
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
#' @importFrom stats coef lm summary.lm
#' @importFrom magrittr %<>% %>%
#' @importFrom dplyr ungroup
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#' net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#' run_zscore(mat, net, minsize=0)
run_zscore <- function(mat,
                    network,
                    .source = source,
                    .target = target,
                    .mor = mor,
                    .likelihood = likelihood,
                    sparse = FALSE,
                    center = FALSE,
                    na.rm = FALSE,
                    minsize = 5L,
                    flavor = "RoKAI"
) {
  
  # NSE vs. R CMD check workaround
  condition <- likelihood <- mor <- p_value <- score <-
    source <- statistic <- target <- NULL
  
  # Check for NAs/Infs in mat
  mat %<>% check_nas_infs
  
  network %>%
    # Convert to standard tibble: source-target-mor.
    rename_net(
      {{ .source }},
      {{ .target }},
      {{ .mor }},
      {{ .likelihood }}
    ) %>%
    filt_minsize(rownames(mat), ., minsize) %>%
    # Preprocessing -------------------------------------------------------
  .fit_preprocessing(mat, center, na.rm, sparse) %>%
    # Model evaluation ----------------------------------------------------
  {.zscore_analysis(.$mat, .$mor_mat, flavor)} %>%
    ungroup()
  
}

#' Wrapper to execute run_zscore() logic on preprocessed data
#'
#' Calculate a z-score from the molecular features of its targets and 
#' normalise it by the background molecular features.
#'
#' @inheritParams run_zscore
#' @param mor_mat
#'
#' @inherit run_zscore return
#' @keywords intern
#' @importFrom stats pnorm
#' @importFrom dplyr filter rename mutate group_by summarise arrange
#' @importFrom tibble rownames_to_column tibble
#' @importFrom tidyr pivot_longer drop_na
#' @importFrom purrr map_dfr
#' @importFrom magrittr %<>% %>%
#' @noRd
.zscore_analysis <- function(mat, mor_mat, flavor) {
  net <- mor_mat %>%
    as.data.frame() %>%
    rownames_to_column("target") %>%
    pivot_longer(cols = -target, values_to = "mor", names_to = "source") %>%
    filter(mor != 0)

  scores <- purrr::map_dfr(seq_len(ncol(mat)), function(exp){
    # Convert column of mat to data frame and drop NA rows
    mat_c <- mat[, exp, drop = FALSE] %>%
      data.frame() %>%
      drop_na() %>%
      rownames_to_column("target")
    
    # Perform inner join and filter NA rows
    KSdata <- full_join(net, mat_c, by = "target") %>%
      drop_na() %>%
      rename("stat" = colnames(.)[4])
    
    # Calculate value based on mor and stat
    KSdata <- KSdata %>%
      mutate(value = mor * stat)
    
    # Aggregate mean values for each source
    Mean.FC <- KSdata %>%
      group_by(source) %>%
      summarise(mS = mean(value), m = n()) %>%
      arrange(source)
    
    # Calculate Enrichment and z-score
    mean_value <- mean(mat_c[, 2], na.rm = TRUE)
    if (flavor == "RoKAI") {
      mean_value <- 0
    }
    Mean.FC <- Mean.FC %>%
      mutate(
        Enrichment = mS / abs(mean(mat_c[,2], na.rm = TRUE)),
        z.score = ((mS - mean_value) * sqrt(m)) / sd(mat_c[, 2], na.rm = TRUE),
        p.value = pnorm(-abs(z.score))
      )
    
    # Prepare and return result as tibble
    tibble(
      statistic = "z_score",
      source = Mean.FC$source,
      condition = colnames(mat_c)[2],
      score = Mean.FC$z.score,
      p_value = Mean.FC$p.value
    )
  })
  # Arrange scores according to the order of sources in mor_mat
  scores <- scores %>%
    arrange(match(source, colnames(mor_mat)))
  
  return(scores)
}

