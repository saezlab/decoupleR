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

# scira and pscira ------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_scira
#'
#' @export
#' @family convert_to_ variants
convert_to_scira <- function(dataset, .source, .target, .target_profile = NULL, clean = FALSE) {
  .check_quos_status({{ .source }}, {{ .target }}, .dots_names = c(".source", ".target"))

  dataset %>%
    transmute_defaults(
      .tf = {{ .source }},
      .target = {{ .target }},
      .mor = {{ .target_profile }},
      .def_col_val = c(.mor = 0)
    ) %>%
    rename(
      tf = .data$.tf,
      target = .data$.target,
      mor = .data$.mor,
    )
}

#' @rdname convert_to_
#'
#' @inheritParams run_pscira
#'
#' @export
#' @family convert_to_ variants
convert_to_pscira <- function(dataset, .source, .target, .target_profile = NULL, clean = FALSE) {
  .check_quos_status({{ .source }}, {{ .target }}, .dots_names = c(".source", ".target"))

  dataset %>%
    transmute_defaults(
      .tf = {{ .source }},
      .target = {{ .target }},
      .mor = {{ .target_profile }},
      .def_col_val = c(.mor = 0)
    ) %>%
    rename(
      tf = .data$.tf,
      target = .data$.target,
      mor = .data$.mor,
    )
}

# mean --------------------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_mean
#'
#' @export
#' @family convert_to_ variants
convert_to_mean <- function(dataset, .source, .target, .mor = NULL, .likelihood = NULL) {
  .check_quos_status({{ .source }}, {{ .target }}, .dots_names = c(".source", ".target"))

  dataset %>%
    transmute_defaults(
      .tf = {{ .source }},
      .target = {{ .target }},
      .mor = {{ .mor }},
      .likelihood = {{ .likelihood }},
      .def_col_val = c(.mor = 0, .likelihood = 1)
    ) %>%
    rename(
      tf = .data$.tf,
      target = .data$.target,
      mor = .data$.mor,
      likelihood = .data$.likelihood
    )
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

#' Stop if any of past quos are missing or NULL.
#'
#' @param ... Quos to evaluate if they are missing or NULL.
#' @param .labels Name corresponding to each quo.
#'
#' @keywords internal
#' @noRd
# TODO be able to use name of dots as name of quo.
.check_quos_status <- function(..., .dots_names) {
  dots <- enquos(...)

  walk2(.x = dots, .y = .dots_names, function(.dot, .name) {
    if (quo_is_missing(.dot)) {
      abort(
        message = str_glue('Quo "{.name}" is missing, with no default.'),
        .subclass = "quo_missing_error"
      )
    }
    if (quo_is_null(.dot)) {
      abort(
        message = str_glue('Quo "{.name}" can not be NULL.'),
        .subclass = "quo_null_error"
      )
    }
  })
}

#' Transmute add column with defaults
#'
#' @description
#' \code{transmute_defaults} combine the \code{\link[dplyr]{transmute}} way of
#' working and combine it with the \code{\link[tibble]{add_column}} to add columns
#' with default values in case they don't exist after transmuting the dataset.
#'
#' @inheritParams dplyr::transmute
#' @param .def_col_val Named vector with columns with default values
#'  if none exist after transmute.
#'
#' @details
#' \code{transmute} adds new variables and drops existing ones.
#' New variables overwrite existing variables of the same name.
#' Variables can be removed by setting their value to NULL.
#'
#' @inherit dplyr::transmute return
#' @keywords internal
#' @import dplyr
#' @import tibble
transmute_defaults <- function(.data, ..., .def_col_val = c()) {
  .data %>%
    transmute(...) %>%
    add_column(., !!!.def_col_val[!names(.def_col_val) %in% names(.)])
}
