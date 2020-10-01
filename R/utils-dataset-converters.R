# main convert_to_ --------------------------------------------------------

#' Convert datasets to standard functions format.
#'
#' @description
#' Convert a dataset to the standard
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

# scira and variants ------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_scira
#'
#' @export
#' @family convert_to_ variants
convert_to_scira <- function(dataset, .source, .target, .profile = NULL, clean = FALSE) {
  .profile <- enquo(.profile)

  if (quo_is_null(.profile)) {
    new_dataset <- dataset %>%
      rename(tf = {{ .source }}, target = {{ .target }}) %>%
      mutate(mor = 0)
  } else {
    new_dataset <- dataset %>%
      rename(tf = {{ .source }}, target = {{ .target }}, mor = {{ .profile }})
  }

  .clean(new_dataset, tf, target, mor, clean = clean)
}

#' @rdname convert_to_
#'
#' @inheritParams run_pscira
#'
#' @export
#' @family convert_to_ variants
convert_to_pscira <- function(dataset, .source, .target, .profile = NULL, clean = FALSE) {
  .profile <- enquo(.profile)

  if (quo_is_null(.profile)) {
    new_dataset <- dataset %>%
      rename(tf = {{ .source }}, target = {{ .target }}) %>%
      mutate(mor = 0)
  } else {
    new_dataset <- dataset %>%
      rename(tf = {{ .source }}, target = {{ .target }}, mor = {{ .profile }})
  }

  .clean(new_dataset, tf, target, mor, clean = clean)
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
