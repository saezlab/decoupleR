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
  get_Rd_metadata <- utils::getFromNamespace (".Rd_get_metadata", "tools")
  dplyr::bind_rows(lapply(db, function(fun){
    name <- get_Rd_metadata(fun, 'name')
    title <- get_Rd_metadata(fun, 'title')
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

#' Pre-processing for methods that fit models
#'
#' - Get only the intersection of target features between `mat` and `network`.
#' - Transform tidy `network` into `matrix` representation with `mor` as value.
#' - If `center` is true, then the expression values are centered by the
#'   mean of expression across the samples.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param sparse Deprecated parameter.
#' @param na.rm Should missing values (including NaN) be omitted from the
#'  calculations of [base::rowMeans()]?
#' @param center Logical value indicating if `mat` must be centered by
#' [base::rowMeans()].
#'
#' @return A named list of matrices to evaluate in methods that fit models, like
#'  `.mlm_analysis()`.
#'  - mat: Features as rows and samples as columns.
#'  - mor_mat: Features as rows and columns as source.
#' @export
.fit_preprocessing <- function(network, mat, center, na.rm, sparse) {
  shared_targets <- intersect(
    rownames(mat),
    network$target
  )
  
  mat <- mat[shared_targets, , drop=FALSE]
  
  mor_mat <- network %>%
    filter(.data$target %in% shared_targets) %>%
    pivot_wider_profile(
      id_cols = .data$target,
      names_from = .data$source,
      values_from = .data$mor,
      values_fill = 0,
      to_matrix = TRUE,
      to_sparse = FALSE
    ) %>%
    .[shared_targets, , drop=FALSE]
  
  likelihood_mat <- network %>%
    filter(.data$target %in% shared_targets) %>%
    pivot_wider_profile(
      id_cols = .data$target,
      names_from = .data$source,
      values_from = .data$likelihood,
      values_fill = 0,
      to_matrix = TRUE
    ) %>%
    .[shared_targets, , drop=FALSE]
  
  weight_mat <- mor_mat * likelihood_mat
  
  if (center) {
    mat <- mat - rowMeans(mat, na.rm)
  }
  
  list(mat = mat, mor_mat = weight_mat)
}
