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
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' filter_regulons(network, .source = tf, min_size = 30, max_size = 50)
filter_regulons <- function(network,
                            .source,
                            min_size = 1,
                            max_size = Inf) {
    network %>%
        add_count({{ .source }}, wt = NULL) %>%
        filter(.data$n >= min_size, .data$n <= max_size) %>%
        select(-.data$n)
}

#' Intersect network target genes with expression matrix.
#'
#' Keep only edges which its target genes belong to the expression matrix.
#'
#' @inheritParams .decoupler_network_format
#' @param .target Maximum number of targets allowed per regulon.
#'
#' @return Filtered tibble.
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' filter_regulons(mat, network, target)
intersect_regulons <- function(mat,
                               network,
                               source,
                               target,
                               minsize
                               ) {
  targets <- rownames(mat)
  network %>%
    filter(!!sym(target) %in% targets) %>%
    group_by(!! sym(source)) %>%
    filter(n() >= minsize)
}
