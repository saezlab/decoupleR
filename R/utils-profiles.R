#' @inherit tidyr::complete title
#'
#' @inherit tidyr::complete description
#'
#' @inheritParams tidyr::complete
#' @inheritParams tidyr::pivot_wider
#' @param sources A named vector or list with the values to expand and get
#'  profile.
#'
#' @return A data frame with the expanded grid of the values passed in
#'  `sources` and filled as specified in the `fill` argument.
#' @examples
#' \dontrun{
#' library(dplyr, warn.conflicts = FALSE)
#' df <- tibble(
#'     group = c(1:2, 1),
#'     item_id = c(1:2, 2),
#'     item_name = c("a", "b", "b"),
#'     value1 = 1:3,
#'     value2 = 4:6
#' )
#'
#' to_get_profile <- list(group = c(1, 2, 3), item_id = c(1, 2))
#'
#' # This will add the combinations of group 3 with the id of the items
#' df %>% get_profile_of(sources = to_get_profile)
#'
#' # You can also choose to fill in missing values
#'
#' # This only fill with "Unknown" the NA values of the column item_name
#' df %>% get_profile_of(
#'     sources = to_get_profile,
#'     values_fill = list(item_name = "Unknown")
#' )
#'
#' # Replace all NAs with "Unkwnon"
#' df %>% get_profile_of(sources = to_get_profile, values_fill = "Unknown")
#' }
#' @keywords internal
#' @seealso [complete][tidyr::complete] [expand][tidyr::expand]
#'
#' @import dplyr
#' @import purrr
#' @import tidyr
get_profile_of <- function(data, sources, values_fill = NA) {
    # The function only allows to reduce or extend the length of the profile,
    # not to add metadata
    stopifnot(all(names(sources) %in% colnames(data)))

    # Drop duplicated entries
    sources <- map(sources, unique)

    # Get combinations of the data and join them to the original data set
    new_data <- lift_dl(expand_grid)(sources) %>%
        left_join(data, by = names(sources))

    if (is_list(values_fill)) {
        replace_na(new_data, replace = values_fill)
    } else if (!is.na(values_fill) && length(values_fill) == 1) {
        new_data %>%
            mutate(across(
                .cols = everything(),
                .fns = ~ replace_na(.x, replace = values_fill)
            ))
    } else {
        new_data
    }
}

#' Pivot a data frame to wider and convert it to matrix
#'
#' @description Generates a kind of table where the rows come from `id_cols`,
#' the columns from `names_from` and the values from `values_from`.
#'
#' @details
#' In the current state of the function, to ensure its operation,
#' the `id_cols` parameter is a single selector.
#'
#' @inheritParams tidyr::pivot_wider
#' @inheritParams tidyr::spread
#' @param to_matrix Logical value indicating if the result should be a matrix.
#'  Parameter is ignored in case `sparse` is `TRUE`.
#' @param to_sparse Logical value indicating whether the resulting matrix
#'  should be sparse or not.
#'
#' @return "widened" data; it is increasing the number of columns and
#'  decreasing the number of rows.
#'
#' @import dplyr
#' @import tibble
#' @import tidyr
#' @importFrom Matrix Matrix
#' @export
#' @examples
#' \dontrun{
#' df <- tibble::tibble(
#'     tf = c("tf_1", "tf_1", "tf_2", "tf_2"),
#'     gene = c("gene_1", "gene_2", "gene_1", "gene_2"),
#'     mor = c(1, -1, 1, -1)
#' )
#'
#' # Return a tibble
#' pivot_wider_profile(
#'     data = df,
#'     id_cols = tf,
#'     names_from = gene,
#'     values_from = mor
#' )
#'
#' # Return a matrix
#' pivot_wider_profile(
#'     data = df,
#'     id_cols = tf,
#'     names_from = gene,
#'     values_from = mor,
#'     to_matrix = TRUE
#' )
#' # Return a sparse Matrix of class "dgCMatrix"
#' pivot_wider_profile(
#'     data = df,
#'     id_cols = tf,
#'     names_from = gene,
#'     values_from = mor,
#'     to_sparse = TRUE
#' )
#' }
#' @keywords internal
pivot_wider_profile <- function(data,
                                id_cols,
                                names_from,
                                values_from,
                                values_fill = NA,
                                to_matrix = FALSE,
                                to_sparse = FALSE,
                                ...) {
    wider_profile <- data %>%
        select({{ id_cols }}, {{ names_from }}, {{ values_from }}) %>%
        pivot_wider(
            id_cols = {{ id_cols }},
            names_from = {{ names_from }},
            values_from = {{ values_from }},
            values_fill = values_fill,
            ...
        ) %>%
        column_to_rownames(var = as_label(enquo(id_cols)))

    if (to_matrix == TRUE || to_sparse == TRUE) {
        if (to_sparse == TRUE) {
            return(Matrix(data = as.matrix(wider_profile), sparse = TRUE))
        } else {
            return(as.matrix(wider_profile))
        }
    }
    wider_profile
}
