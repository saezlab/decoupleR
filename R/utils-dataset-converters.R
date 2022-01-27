#' Rename network
#'
#' @description
#' Renames a given network to these column names: .source, .target, .mor, If 
#' .mor is not provided, then the function sets them to default values.
#' 
#' @inheritParams .decoupler_network_format
#' @param def_mor Default value for .mor when not provided.
#' 
#' @export
#' @examples 
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' rename_net(network, tf, target, mor)
rename_net <- function(network,
                       .source,
                       .target,
                       .mor = NULL,
                       .likelihood = NULL,
                       def_mor = 1) {
    
    .check_quos_status({{ .source }}, {{ .target }}, 
                       .dots_names = c(".source", ".target"))
    if (!'likelihood' %in% colnames(network)){
        network <- network %>% mutate(likelihood=1)
    }
    network <- network %>%
        convert_f_defaults(
            source = {{ .source }},
            target = {{ .target }},
            mor = {{ .mor }},
            likelihood = {{ .likelihood }},
            .def_col_val = c(mor = def_mor, likelihood=1)
        )
    if (any(network$likelihood != 1)) {
        warning(".likelihood argument is deprecated, it will be set to 1. From now
                on, weights of regulation should go into the .mor column.")
    }
    check_repeated_edges(network)
    network <- network %>% mutate(likelihood=1)
    network
}

#' Extract sets
#'
#' @description
#' Extracts feature sets from a renamed network (see [decoupleR::rename_net]).
#' 
#' @inheritParams .decoupler_network_format
#' 
#' @export
#' 
#' @examples 
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' network <- rename_net(network, tf, target, mor, likelihood)
#' extract_sets(network)
extract_sets <- function(network) {
    network %>%
        group_by(.data$source) %>%
        summarise(
            regulons = set_names(list(.data$target), .data$source[1]),
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
            rlang::abort(
                message = stringr::str_glue(
                    'Quo "{.name}" is missing, with no default.'
                ),
                class = "quo_missing_error"
            )
        }
        if (quo_is_null(.dot)) {
            rlang::abort(
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

        rlang::abort(
            message = stringr::str_glue(
                "Output columns are different than expected.\n",
                "Expected: {expected_columns}\n",
                "Extra: {extra_cols}\n",
                "Removed: {removed_cols}"
            ),
            class = "different_set_columns"
        )
    }

    .data
}

#' Check if network contains repeated edges
#'
#' @param network Network in tibble format.
#' @noRd
check_repeated_edges <- function(network){
    repeated <- network %>%
        group_by(.data$source, .data$target) %>%
        filter(n()>1)
    if (nrow(repeated) > 1){
        stop('Network contains repeated edges, please remove them.')
    }
}

#' Check if mat contains Nans or Infs
#'
#' @param mat Matrix in matrix format.
#' @noRd
check_nas_infs <- function(mat){
    mat <- as.matrix(mat)
    if (any(is.infinite(mat) | is.na(mat))){
        stop('Mat contains NAs or Infs, please remove them.')
    }
}
