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
#' @param minsize Integer indicating the minimum number of targets per source.
#'
#' @return Filtered network.
#' @export
#' @examples 
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' network <- rename_net(network, tf, target, mor, likelihood)
#' filt_minsize(rownames(mat), network, minsize = 5)
filt_minsize <- function(mat_f_names, network, minsize = 5){
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
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#' network <- rename_net(network, tf, target, mor, likelihood)
#' .fit_preprocessing(network, mat, center = FALSE, na.rm = FALSE, sparse = FALSE)
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
  
  if (is.null(colnames(mat))){
    colnames(mat) <- 1:ncol(mat)
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
                        .likelihood = NULL){
  
  source <- as.symbol(.source)
  target <- as.symbol(.target)
  mor <- as.symbol(.mor)
  
  network <- network %>%
    dplyr::mutate(likelihood=1) %>%
    dplyr::mutate(weight = (!!mor)*.data$likelihood) %>% 
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

#' Generate a toy `mat` and `network`.
#' 
#' @param n_samples Number of samples to simulate.
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
#'
#' @return List containing `mat` and `network`.
#' @export
#' @examples
#' data <- get_toy_data()
#' mat <- data$mat
#' network <- data$network
get_toy_data <- function(n_samples = 24, seed = 42){
  network <- tibble::tibble(
    source = c('T1','T1','T1','T2','T2','T2','T3','T3','T3','T3'),
    target = c('G01','G02','G03','G06','G07','G08','G06','G07','G08','G11'),
    mor = c(1,1,0.7,1,0.5,1,-0.5,-3,-1,1)
  )
  
  n_features <- 12
  n <- round(n_samples/2)
  res = n_samples %% 2
  r1 <- c(8,8,8,8,0,0,0,0,0,0,0,0)
  r2 <- c(0,0,0,0,8,8,8,8,0,0,0,0)
  rep(matrix(r1, ncol=12),n)
  matrix(c(rep(r1, n), ncol=12, nrow=12))
  matrix(c(rep(r1, n),rep(r2, n)), ncol = 24)
  mat <- matrix(c(rep(r1,n),rep(r2,n)), nrow=12)
  withr::with_seed(seed, {
    rand <- matrix(abs(stats::rnorm(dim(mat)[1]*dim(mat)[2])), nrow=12)
  })
  mat <- mat + rand
  rownames(mat) <- c('G01','G02','G03','G04','G05','G06','G07','G08','G09','G10','G11','G12')
  colnames(mat) <- lapply(1:dim(mat)[2], function(i){paste0('S',sprintf("%02d", i))})
  
  return(list(mat=mat, network=network))
}
