# main convert_to_ --------------------------------------------------------

#' Convert a network to run under the method of interest.
#'
#' @description
#' Convert a long-format network to the suggested standard for the
#' specified `run_{statistic}()`. If the default parameters are not modified,
#' then the function sets its own null values for those columns.
#'
#' @inheritParams .decoupler_network_format
#'
#' @return
#'
#' + `convert_to_`
#'    Return same as input.
#' * `convert_to_gsva()`
#'    Return a list of regulons suitable for [GSVA::gsva()].
#' * `convert_to_mean()`
#'    Return a tibble with four columns: `tf`, `target`, `mor` and `likelihood`.
#' * `convert_to_ora()`
#'    Return a named list of regulons; tf with associated targets.
#' * `convert_to_pscira()`
#'    Returns a tibble with three columns: `tf`, `target` and `mor`.
#' * `convert_to_scira()`
#'    Returns a tibble with three columns: `tf`, `target` and `mor`.
#' * `convert_to_viper()`
#'    Return a list of regulons suitable for [viper::viper()]
#'
#' @name convert_to_
#' @rdname convert_to_
#' @family convert_to_ variants
#'
#' @seealso [convert_f_defaults()]
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' convert_to_(network)
#' convert_to_gsva(network, tf, target)
#' convert_to_mean(network, tf, target, mor, likelihood)
#' convert_to_ora(network, tf, target)
#' convert_to_pscira(network, tf, target, mor)
#' convert_to_scira(network, tf, target, mor)
#' convert_to_viper(network, tf, target, mor, likelihood)
convert_to_ <- function(network) invisible(network)

# scira and pscira ------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_scira
#'
#' @family convert_to_ variants
#' @export
convert_to_scira <- function(network, .source, .target, .mor = NULL) {
    .check_quos_status({{ .source }}, {{ .target }},
        .dots_names = c(".source", ".target")
    )

    network %>%
        convert_f_defaults(
            tf = {{ .source }},
            target = {{ .target }},
            mor = {{ .mor }},
            .def_col_val = c(mor = 0)
        ) %>%
        mutate(mor = sign(.data$mor))
}

#' @rdname convert_to_
#'
#' @inheritParams run_pscira
#'
#' @family convert_to_ variants
#' @export
convert_to_pscira <- function(network, .source, .target, .mor = NULL) {
    .check_quos_status({{ .source }}, {{ .target }},
        .dots_names = c(".source", ".target")
    )

    network %>%
        convert_f_defaults(
            tf = {{ .source }},
            target = {{ .target }},
            mor = {{ .mor }},
            .def_col_val = c(mor = 0)
        ) %>%
        mutate(mor = sign(.data$mor))
}

# mean --------------------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_mean
#'
#' @family convert_to_ variants
#' @export
convert_to_mean <- function(network,
                            .source,
                            .target,
                            .mor = NULL,
                            .likelihood = NULL) {
    .check_quos_status({{ .source }}, {{ .target }},
        .dots_names = c(".source", ".target")
    )

    network %>%
        convert_f_defaults(
            tf = {{ .source }},
            target = {{ .target }},
            mor = {{ .mor }},
            likelihood = {{ .likelihood }},
            .def_col_val = c(mor = 0, likelihood = 1)
        ) %>%
        mutate(mor = sign(.data$mor))
}

# viper -------------------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_viper
#'
#' @family convert_to_ variants
#' @export
convert_to_viper <- function(network,
                             .source,
                             .target,
                             .mor = NULL,
                             .likelihood = NULL) {
    .check_quos_status({{ .source }}, {{ .target }},
        .dots_names = c(".source", ".target")
    )

    network %>%
        convert_f_defaults(
            tf = {{ .source }},
            target = {{ .target }},
            mor = {{ .mor }},
            likelihood = {{ .likelihood }},
            .def_col_val = c(mor = 0, likelihood = 1)
        ) %>%
        mutate(mor = sign(.data$mor)) %>%
        split(.$tf) %>%
        map(~ {
            list(
                tfmode = set_names(.x$mor, .x$target),
                likelihood = .x$likelihood
            )
        })
}

# gsva --------------------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_gsva
#'
#' @family convert_to_ variants
#' @export
convert_to_gsva <- function(network, .source, .target) {
    .check_quos_status({{ .source }}, {{ .target }},
        .dots_names = c(".source", ".target")
    )

    network %>%
        convert_f_defaults(
            tf = {{ .source }},
            target = {{ .target }}
        ) %>%
        group_by(.data$tf) %>%
        summarise(
            regulons = set_names(list(.data$target), .data$tf[1]),
            .groups = "drop"
        ) %>%
        pull(.data$regulons)
}

# ora ---------------------------------------------------------------------

#' @rdname convert_to_
#'
#' @inheritParams run_ora
#'
#' @family convert_to_ variants
#' @export
convert_to_ora <- function(network, .source, .target) {
    .check_quos_status({{ .source }}, {{ .target }},
        .dots_names = c(".source", ".target")
    )

    network %>%
        convert_f_defaults(
            tf = {{ .source }},
            target = {{ .target }}
        ) %>%
        group_by(.data$tf) %>%
        summarise(
            regulons = set_names(list(.data$target), .data$tf[1]),
            .groups = "drop"
        ) %>%
        pull(.data$regulons)
}

# Helper functions --------------------------------------------------------

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
                message = stringr::str_glue(
                    'Quo "{.name}" is missing, with no default.'
                ),
                class = "quo_missing_error"
            )
        }
        if (quo_is_null(.dot)) {
            abort(
                message = stringr::str_glue('Quo "{.name}" can not be NULL.'),
                class = "quo_null_error"
            )
        }
    })
}

#' Rename columns and add defaults values if column not present
#'
#' @description
#' `convert_f_defaults()` combine the [dplyr::rename()] way of
#' working and with the [tibble::add_column()] to add columns
#' with default values in case they don't exist after renaming data.
#'
#' @inheritParams dplyr::rename
#' @param .def_col_val Named vector with columns with default values
#'  if none exist after rename.
#' @param .use_dots Should a dot prefix be added to renamed variables?
#' This will allow swapping of columns.
#'
#' @details
#' The objective of using .use_dots is to be able to swap columns which,
#' by default, is not allowed by the [dplyr::rename()] function.
#' The same behavior can be replicated by simply using the [dplyr::select()],
#' however, the select evaluation allows much more flexibility so that
#' unexpected results could be obtained. Despite this, a future implementation
#' will consider this form of execution to allow renaming the same
#' column to multiple ones (i.e. extend dataframe extension).
#'
#' @return
#' An object of the same type as .data. The output has the following properties:
#' - Rows are not affected.
#' - Column names are changed.
#' - Column order is the same as that of the function call.
#' @export
#' @importFrom tidyselect eval_rename
#' @examples
#'
#' df <- tibble::tibble(x = 1, y = 2, z = 3)
#'
#' # Rename columns
#' df <- tibble::tibble(x = 1, y = 2)
#' convert_f_defaults(
#'     .data = df,
#'     new_x = x,
#'     new_y = y,
#'     new_z = NULL,
#'     .def_col_val = c(new_z = 3)
#' )
convert_f_defaults <- function(.data,
                               ...,
                               .def_col_val = c(),
                               .use_dots = TRUE) {
    expected_columns <- match.call(expand.dots = FALSE)$... %>%
        names() %>%
        unique()

    .expr <- expr(c(...))
    if (.use_dots) .expr <- expr(c(. = !!.expr))

    # Return rename changes with dot prefix variables.
    loc <- eval_rename(.expr, data = .data)

    .data %>%
        select(all_of(loc)) %>%
        {
            # Remove prefix dots generated by eval_rename()
            if (.use_dots) {
                rename_with(., ~ stringr::str_remove(.x, "...."))
            } else {
                .
            }
        } %>%
        add_column(., !!!.def_col_val[!names(.def_col_val) %in% names(.)]) %>%
        .check_expected_columns(expected_columns = expected_columns)
}

#' Check if data contains specific columns
#'
#' If `.data` present more or less columns than expected
#' then the function will abort execution, otherwise it will
#' return the same input data.
#'
#' @inheritParams convert_f_defaults
#' @param expected_columns Name of the columns that must make a total match
#'  with the expected columns
#'
#' @return `.data`
#'
#' @noRd
.check_expected_columns <- function(.data, expected_columns) {
    # Get data columns.
    data_cols <- names(.data)

    # Calculate symmetric difference
    diff_cols <- setdiff(
        x = union(expected_columns, data_cols),
        y = intersect(expected_columns, data_cols)
    )

    # Abort execution if there is an inconsistency in the output results
    if (!is_empty(diff_cols)) {
        extra_cols <- setdiff(diff_cols, expected_columns) %>%
            paste(collapse = ", ")
        removed_cols <- intersect(expected_columns, diff_cols) %>%
            paste(collapse = ", ")
        expected_columns <- paste(expected_columns, collapse = ", ")

        out_message <- stringr::str_glue(
            "Output columns are different than expected.\n",
            "Expected: {expected_columns}\n",
            "Extra: {extra_cols}\n",
            "Removed: {removed_cols}"
        )

        abort(message = out_message, class = "different_set_columns")
    }

    .data
}
