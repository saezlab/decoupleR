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

#' Filter sources with minsize targets
#'
#' Filter sources of a net with less than minsize targets
#'
#' @param mat_f_names Feature names of mat.
#' @inheritParams .decoupler_network_format
#'
#' @return Filtered network.
#' @export
filt_minsize <- function(mat_f_names, network, minsize=5){
  # Find shared targets
  shared_targets <- intersect(
    mat_f_names,
    network$target
  )
  
  # Find sizes of sources after intersect and filter by minsize
  sources <- network %>%
    dplyr::filter(.data$target %in% shared_targets) %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(n=dplyr::n()) %>%
    dplyr::filter(.data$n >= minsize) %>%
    dplyr::pull(.data$source)
  
  # Filter sources
   network <- network %>%
    dplyr::filter(.data$source %in% sources)
   
   if (nrow(network) == 0) {
     stop(stringr::str_glue('Network is empty after intersecting it with mat and
     filtering it by sources with at least {minsize} targets. Make sure mat and 
     network have shared target features or reduce the number assigned to minsize'))
   }
   return(network)
}

#' Pre-processing for methods that fit networks.
#'
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
  # Create empty mor_mat from original feature universe from mat, then fill in
  sources <- unique(network$source)
  targets <- rownames(mat)
  mor_mat <- matrix(0, ncol = length(sources), nrow=nrow(mat))
  colnames(mor_mat) <- sources
  rownames(mor_mat) <- targets
  weights <- network$mor * network$likelihood
  for (i in 1:nrow(network)) {
    .source <- network$source[[i]]
    .target <- network$target[[i]]
    .weight <- weights[[i]]
    if (.target %in% targets) {
      mor_mat[[.target,.source]] <- .weight
    }
  }
  
  if (center) {
    mat <- mat - rowMeans(mat, na.rm)
  }
  
  list(mat = mat, mor_mat = mor_mat)
}

#' Check correlation (colinearity) 
#'
#' Checks the correlation across the regulators in a network.
#' @inheritParams .decoupler_network_format
#'
#' @return Correlation pairs tibble.
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' check_corr(network, .source='tf')
check_corr <- function(network, 
                        .source = "source",
                        .target = "target", 
                        .mor = "mor", 
                        .likelihood = "likelihood"){
  
  source <- as.symbol(.source)
  target <- as.symbol(.target)
  mor <- as.symbol(.mor)
  likelihood <- as.symbol(.likelihood)
  
  network <- network %>% 
    dplyr::mutate(weight = (!!mor)*(!!likelihood)) %>% 
    dplyr::select(!!source, !!target, .data$weight)
  network_wide <- network %>%
    tidyr::pivot_wider(names_from = !!target, values_from = .data$weight, values_fill = 0) %>%
    tibble::column_to_rownames(.source)
  
  cor_source <- stats::cor(t(network_wide))
  cor_source[lower.tri(cor_source, diag = TRUE)] <- NA
  
  cor_source <- cor_source %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(.source) %>%
    tidyr::pivot_longer(!(!!source), names_to = paste0(.source, ".2"), values_to = "correlation") %>%
    dplyr::filter(!is.na(.data$correlation)) %>% 
    dplyr::arrange(desc(.data$correlation))
  cor_source
}

