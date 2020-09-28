#' @inherit tidyr::complete title
#'
#' @inherit tidyr::complete description
#'
#' @inheritParams tidyr::complete
#' @inheritParams tidyr::pivot_wider
#' @param sources A named vector or list with the values to expand and get profile.
#'
#' @return A data frame with the expanded grid of the values passed in
#'  \code{sources} and filled as specified in the \code{fill} argument.
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' df <- tibble(
#'   group = c(1:2, 1),
#'   item_id = c(1:2, 2),
#'   item_name = c("a", "b", "b"),
#'   value1 = 1:3,
#'   value2 = 4:6
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
#' df %>% get_profile_of(sources = to_get_profile, values_fill = list(item_name = "Unknown"))
#'
#' # Replace all NAs with "Unkwnon"
#' df %>% get_profile_of(sources = to_get_profile, values_fill = "Unknown")
#' @seealso \link[tidyr]{complete} \link[tidyr]{expand}
#'
#' @importFrom tidyr complete expand_grid replace_na replace_na
#' @importFrom dplyr left_join mutate across everything
#' @importFrom purrr map lift_dl is_list
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
      mutate(across(everything(), ~ replace_na(.x, replace = values_fill)))
  } else {
    new_data
  }
}
