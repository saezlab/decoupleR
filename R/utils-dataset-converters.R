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
    convert_f_defaults(
      tf = {{ .source }},
      target = {{ .target }},
      mor = {{ .target_profile }},
      .def_col_val = c(mor = 0)
    ) %>%
    mutate(mor = sign(.data$mor))
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
    convert_f_defaults(
      tf = {{ .source }},
      target = {{ .target }},
      mor = {{ .target_profile }},
      .def_col_val = c(mor = 0)
    ) %>%
    mutate(mor = sign(.data$mor))
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
    convert_f_defaults(
      tf = {{ .source }},
      target = {{ .target }},
      mor = {{ .mor }},
      likelihood = {{ .likelihood }},
      .def_col_val = c(mor = 0, likelihood = 1)
    ) %>%
    mutate(mor = sign(.data$mor))
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

#' Rename and add column with defaults
#'
#' @description
#' \code{rename_defaults} combine the \code{\link[dplyr]{rename}} way of
#' working and with the \code{\link[tibble]{add_column}} to add columns
#' with default values in case they don't exist after renaming the dataset.
#'
#' @inheritParams dplyr::rename
#' @param .def_col_val Named vector with columns with default values
#'  if none exist after rename.
#'
#' @inherit dplyr::rename return
#' @keywords internal
#' @import dplyr
#' @import tibble
rename_defaults <- function(.data, ..., .def_col_val = c()) {
  .data %>%
    rename(...) %>%
    add_column(., !!!.def_col_val[!names(.def_col_val) %in% names(.)])
}


#' Rename columns and add defaults values if column not present
#'
#' @description
#' \code{rename_defaults} combine the \code{\link[dplyr]{rename}} way of
#' working and with the \code{\link[tibble]{add_column}} to add columns
#' with default values in case they don't exist after renaming the dataset.
#'
#' @inheritParams dplyr::rename
#' @param .def_col_val Named vector with columns with default values
#'  if none exist after rename.
#' @param .use_dots Should a dot prefix be added to renamed variables?
#' This will allow swapping of columns.
#'
#' @details
#' The objective of using .use_dots is to be able to swap columns which, by default,
#' is not allowed by the "rename" function. The same behavior can be replicated
#' by simply using the "select" function, however, the select evaluation allows
#' much more flexibility so that unexpected results could be obtained.
#' Despite this, a future implementation will consider this form of execution
#' to allow renaming the same column to multiple ones (i.e. extend dataframe extension).
#'
#' @return An object of the same type as .data. The output has the following properties:
#' \itemize{
#'   \item Rows are not affected.
#'   \item Column names are changed.
#'   \item Column order is the same as that of the function call.
#' }
#' @export
#' @importFrom tidyselect eval_rename
convert_f_defaults <- function(.data,
                               ...,
                               .def_col_val = c(),
                               .use_dots = TRUE) {
  out_cols <- match.call(expand.dots = FALSE)$... %>%
    names() %>%
    unique()

  .expr <- expr(c(...))
  if (.use_dots) .expr <- expr(c(. = !!.expr))

  # Return rename changes with dot prefix variables.
  loc <- eval_rename(expr(c(. = c(...))), data = .data)

  .data <- .data %>%
    select(loc) %>%
    {
      # Remove prefix dots generated by eval_rename()
      if (.use_dots) {
        rename_with(., ~ str_remove(.x, "...."))
      } else {
        .
      }
    } %>%
    add_column(., !!!.def_col_val[!names(.def_col_val) %in% names(.)])

  # Check output data columns
  data_cols <- names(.data)

  # Calculate symmetric difference
  diff_cols <- setdiff(
    x = union(out_cols, data_cols),
    y = intersect(out_cols, data_cols)
  )

  # Abort execution if there is an inconsistency in the output results
  if (!is_empty(diff_cols)) {
    extra_cols <- setdiff(diff_cols, out_cols) %>%
      paste(collapse = ", ")
    removed_cols <- intersect(out_cols, diff_cols) %>%
      paste(collapse = ", ")
    out_cols <- paste(out_cols, collapse = ", ")

    out_message <- str_glue(
      "Output columns are different than expected.\n",
      "Expected: {out_cols}\n",
      "Extra: {extra_cols}\n",
      "Removed: {removed_cols}"
    )

    abort(message = out_message, class = "assasassas")
  }

  .data
}
