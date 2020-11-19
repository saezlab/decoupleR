#' VIPER wrapper
#'
#' This function is a convenient wrapper for the [viper::viper()] function.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param options A list of named options to pass to
#' [viper::viper()] such as `minsize` or `method`.
#'  These options should not `include`, `eset` or `regulon`.
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the conditions. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#'  5. `statistic_time`: Internal execution time indicator.
#' @export
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import purrr
#' @import tidyr
#' @import viper
run_viper <- function(mat,
                      network,
                      .source = .data$tf,
                      .target = .data$target,
                      .mor = .data$mor,
                      .likelihood = .data$likelihood,
                      options = list()) {

  # Before to start ---------------------------------------------------------
  .start_time <- Sys.time()

  network <- network %>%
    convert_to_viper({{ .source }}, {{ .target }}, {{ .mor }}, {{ .likelihood }})

  # Analysis ----------------------------------------------------------------
  exec(
    .fn = viper::viper,
    eset = mat,
    regulon = network,
    !!!options
  ) %>%
    as.data.frame() %>%
    rownames_to_column("tf") %>%
    pivot_longer(-.data$tf, names_to = "condition", values_to = "score") %>%
    add_column(statistic = "viper", .before = 1) %>%
    mutate(statistic_time = difftime(Sys.time(), .start_time))
}
