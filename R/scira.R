#' SCIRA (Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculates TF activity according to \href{https://www.biorxiv.org/content/10.1101/553040v1.full.pdf}{SCIRA}.
#'
#' @details
#' Estimation of regulatory activity: A linear regression of the expression profile
#' is performed against the "target profile" of the given TF, where in the target
#' profile, any member of the regulon will be assigned a +1 for activated actions,
#' a -1 for inhibitory activations and, finally, a 0 for all genes that are not
#' a member of the TF regulon. TF activity is then defined as the t-statistic of
#' this linear regression.
#'
#' @param emat A named expression matrix. E.g Genes in rows and samples in columns.
#' @param regulons A data frame of regulons in table format.
#' @param .source Column name in regulons with the factors to calculate the enrichment scores.
#' @param .target Column name in regulons that relates the labels of the rows in the expression matrix.
#' @param .profile Column name in regulons that contains the target profile.
#' @param .sparse Logical value indicating if the generated profile matrix should be sparse.
#'
#' @return A long format tibble of the enrichment results for each set of genes
#'  across the samples. Resulting tibble contains the following columns:
#'  \enumerate{
#'    \item \code{source} Resources from the \code{.source} column of the \code{regulons} data frame.
#'    \item \code{condition} Samples representing each column of \code{emat}.
#'    \item \code{score} Regulatory activity of each resource for each condition.
#'  }
#' @export
#'
#' @importFrom tibble tibble
#' @importFrom stats lm summary.lm coef
#' @importFrom dplyr pull rename rowwise mutate
#' @importFrom purrr lift_dl
#' @importFrom tidyr expand_grid
run_scira <- function(emat,
                      regulons,
                      .source = .data$tf,
                      .target = .data$target,
                      .profile = .data$mor,
                      .sparse = FALSE) {

  # Internal functions ------------------------------------------------------

  # Map model data
  #
  # Build a data set with the necessary values to evaluate the model.
  #
  # @param source Tf index to extract the corresponding associated profiles.
  # @param condition Sample index to extract its expression values.
  #
  # @return Tibble with the data to use to evaluate the model.
  map_model_data <- function(source, condition) {
    tibble(expression = emat[, condition], profile = profile_mat[source, ])
  }

  # Evaluate model
  #
  # Calculate the regulatory activity of a tf under a given condition.
  #
  # @param data A data frame, data frame extension (e.g. a tibble).
  #
  # @return t-value corresponding to beta 1 parameter of linear regression.
  evaluate_model <- function(data) {
    t_values <- lm(expression ~ profile, data = data) %>%
      summary.lm() %>%
      coef() %>%
      .[, "t value"]

    out <- tryCatch(
      {
        return(t_values[[2]]) # Using the t value of the beta1 coefficient.
      },
      error = function(error) {
        # message("t-value value of the coefficient B1 of the linear regression was not accessible.")
        return(NA)
      }
    )
    return(out)
  }

  # Preprocessing -----------------------------------------------------------

  # Extract labels that will map to the expression and profile matrices
  sources <- regulons %>%
    pull({{ .source }}) %>%
    unique()

  conditions <- colnames(emat)

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  profile_mat <- regulons %>%
    rename(source = {{ .source }}, target = {{ .target }}, profile = {{ .profile }}) %>%
    get_profile_of(
      sources = list(source = sources, target = rownames(emat)),
      values_fill = list(profile = 0)
    ) %>%
    pivot_wider_profile(.data$source, .data$target, .data$profile, to_matrix = TRUE, to_sparse = .sparse)


  # Model evaluation --------------------------------------------------------

  # Allocate the space for all combinations of sources and conditions
  # and evaluate the proposed model.
  lift_dl(expand_grid)(list(source = sources, condition = conditions)) %>%
    rowwise() %>%
    mutate(score = map_model_data(.data$source, .data$condition) %>%
      evaluate_model())
}

#' PSCIRA (Permutation Single Cell Inference of Regulatory Activity)
#'
#' @description
#' Calculates TF activity according to \href{https://www.biorxiv.org/content/10.1101/553040v1.full.pdf}{SCIRA}.
#'
#' @inherit run_scira details
#'
#' @inheritParams run_scira
#' @param times Number of replications.
#' @param seed A single value, interpreted as an integer, or NULL.
#'
#' @inherit run_scira return
#' @export
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr pull rename group_by summarise left_join mutate transmute
#' @importFrom purrr map_dfr
run_pscira <- function(emat,
                       regulons,
                       .source = .data$tf,
                       .target = .data$target,
                       .profile = .data$mor,
                       .sparse = TRUE,
                       times = 10,
                       seed = 42) {

  # Internal functions ------------------------------------------------------

  # Map model data
  #
  # Build a data set with the necessary values to evaluate the model.
  #
  # @param .emat Expression matrix.
  # @param random Logical value that indicates whether the rows of the matrix should be shuffled or not.
  #
  # @return Expression matrix.
  map_model_data <- function(.emat, random = FALSE) {
    if (random) {
      return(.emat[sample(nrow(.emat)), ])
    } else {
      .emat
    }
  }

  # Evaluate model
  #
  # Calculates the regulatory activity of all tfs with respect to its
  # associated profile for each condition.
  #
  # @param .emat Expression matrix.
  #
  # @return Tibble with tf regulatory activity for each tf-sample pair.
  evaluate_model <- function(.emat) {
    (profile_mat %*% .emat) %>%
      as.matrix() %>%
      as.data.frame() %>%
      rownames_to_column("source") %>%
      pivot_longer(-.data$source, names_to = "condition")
  }

  # Preprocessing -----------------------------------------------------------

  # Extract labels that will map to the expression and profile matrices
  sources <- regulons %>%
    pull({{ .source }}) %>%
    unique()

  # Ensures column matching, expands the target profile to encompass all targets
  # in the expression matrix for each source, and converts the result to a matrix.
  profile_mat <- regulons %>%
    rename(source = {{ .source }}, target = {{ .target }}, profile = {{ .profile }}) %>%
    get_profile_of(
      sources = list(source = sources, target = rownames(emat)),
      values_fill = list(profile = 0)
    ) %>%
    pivot_wider_profile(.data$source, .data$target, .data$profile, to_matrix = TRUE, to_sparse = .sparse)

  # Convert to matrix to ensure that matrix multiplication works
  # in case emat is a labelled dataframe.
  emat <- as.matrix(emat)

  # Evaluate model ----------------------------------------------------------

  set.seed(seed)
  map_dfr(1:times, ~ map_model_data(emat, random = TRUE) %>%
    evaluate_model()) %>%
    group_by(.data$source, .data$condition) %>%
    summarise(.mean = mean(.data$value), .sd = stats::sd(.data$value), .groups = "drop") %>%
    left_join(evaluate_model(emat), by = c("source", "condition")) %>%
    mutate(score = (.data$value - .data$.mean) / .data$.sd) %>%
    transmute(.data$source, .data$condition, .data$score)
}
