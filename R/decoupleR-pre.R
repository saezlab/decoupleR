#' Filter network by size of regulons
#'
#' Keep only sources which satisfied the condition `min_size >= n <= max_size`,
#' where `n` denotes the number of targets per source.
#'
#' @inheritParams .decoupler_network_format
#' @param min_size Minimum number of targets allowed per regulon.
#' @param max_size Maximum number of targets allowed per regulon.
#'
#' @return Filtered tibble.
#' @export
filter_regulons <- function(
    network,
    .source,
    min_size = 1,
    max_size = Inf) {

    network %>%
        add_count({{ .source }}, wt = NULL) %>%
        filter(.data$n >= min_size, .data$n <= max_size) %>%
        select(-.data$n)
}
