#' Show methods
#'
#' Prints the methods available in decoupleR. The first column correspond to
#' the function name in decoupleR and the second to the method's full name.
#'
#' @export
#' @examples
#' show_methods()
show_methods <- function(){
  db <- tools::Rd_db("decoupleR")
  db <- db[grep("run_*", names(db), value = TRUE)]
  dplyr::bind_rows(lapply(db, function(fun){
    name <- tools:::.Rd_get_metadata(fun, 'name')
    title <- tools:::.Rd_get_metadata(fun, 'title')
    tibble::tibble(Function=name, Name=title)
  }))
}

#' Intersect network target features with input matrix.
#'
#' Keep only edges which its target features belong to the input matrix.
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param minsize Minimum number of targets per source allowed.
#'
#' @return Filtered tibble.
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' intersect_regulons(mat, network, tf, target, minsize=5)
intersect_regulons <- function(mat,
                               network,
                               .source,
                               .target,
                               minsize
) {
  .source<- as.name(substitute(.source))
  .target<- as.name(substitute(.target))
  .source <- enquo(.source)
  .target <- enquo(.target)
  targets <- rownames(mat)
  network %>%
    filter(!!.target %in% targets) %>%
    group_by(!!.source) %>%
    filter(n() >= minsize)
}
