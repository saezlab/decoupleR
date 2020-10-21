# main convert_to_ --------------------------------------------------------

#' Convert data sets to run under the method of interest.
#'
#' @description
#' Convert the data set to the suggested standard for the specified function.
#' If the default parameters are not modified, then the function sets its own
#' null values for those columns.
#'
#' @param dataset A data frame or data frame extension (e.g. a tibble) to convert.
#' @param clean Logical value that indicates whether to keep only the columns
#'  necessary to run the function \code{convert_to_[fun]}.
#'
#' @return Returns a tibble with the necessary columns to evaluate the method
#'  to which the dataset is being converted.
#'
#' @name convert_to_
#' @rdname convert_to_
#' @family convert_to_ variants
#'
#' @export
convert_to_ <- function(dataset, clean) invisible()



#' @rdname convert_to_
#'
#' @inheritParams run_singscore
#'
#' @export
#' @family convert_to_ variants
convert_to_singscore <- function(dataset,
                                 .source,
                                 .target,
                                 .target_profile = NULL,
                                 minsize) {

  .missing_quos({{ .source }}, {{ .target }}, .labels = c(".source", ".target"))

  .target_profile <- enquo(.target_profile)

  if (quo_is_null(.target_profile)) {
    new_dataset <- dataset %>%
      dplyr::mutate(mor = 0) %>%
      dplyr::mutate(likelihood = 0) %>%
      dplyr::rename(geneset = {{ .source }}, gene = {{ .target }}) %>%
      as_tibble() %>%
      undirected2singscore(., minsize)
  } else {
    new_dataset <- dataset %>%
      dplyr::rename(geneset = {{ .source }},
                    gene = {{ .target }},
                    mor = {{ .target_profile }}) %>%
      mutate(mor=sign(mor)) %>%
      directed2singscore(., minsize)
  }
}


# Helper functions --------------------------------------------------------

#' Remove unnecessary variables
#'
#' @inheritParams convert_to_
#' @inheritParams dplyr::select
#'
#' @keywords internal
#' @noRd
.clean <- function(dataset, ..., clean) {
  if (clean) {
    select(dataset, ...)
  } else {
    dataset
  }
}

#' Stop if argument in function is missing.
#'
#' @param ... a
#' @param .labels a
#'
#' @keywords internal
#' @noRd
.missing_quos <- function(..., .labels) {
  vars <- enquos(...)
  walk2(.x = vars, .y = .labels, function(.var, .label) {
    if (quo_is_missing(.var)) {
      stop(str_glue('argument "{.label}" is missing, with no default'))
    }
  })
}
