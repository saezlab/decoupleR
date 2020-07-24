#' VIPER wrapper
#'
#' This function is a convenient wrapper for the
#' \code{\link[=viper]{viper::viper()}} function.
#'
#' @param emat An expression matrix with genes (HGNC symbol) in rows and samples
#'  in columns.
#' @param genesets A data frame of gene sets. The structure is dependent on the
#' gene set resource.
#' @param options A list of named options to pass to
#' \code{\link[=viper]{viper::viper()}} such as \code{minsize} or \code{method}.
#'  These options should not \code{include}, \code{eset} or \code{regulon}.
#' @param gs_resource A character indicating the gene set resource. Based on
#' this argument \code{genesets} is processed futher via diverse helper
#' functions (e.g. \code{\link[=dorothea2viper]{dorothea2viper()}}).
#' @param tidy Logical, whether computed viper scores should be returned in a
#' tidy format.
#'
#' @return A matrix of normalized enrichment scores for each gene set across all
#'  samples. If \code{tidy} is TRUE the normalized enrichment scores are retured
#'  in a tidy format.
#'
#' @export
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import purrr
#' @import tidyr
#' @import viper
run_viper = function(emat, genesets, options = list(), gs_resource, tidy = F) {

  genesets = switch(gs_resource,
    "dorothea" = dorothea2viper(genesets),
    "progeny" = progeny2viper(genesets))

  viper_res = do.call(viper,
                      c(list(eset = emat,
                             regulon = make_viper_genesets(genesets)),
                        options))

  if (tidy) {
    metadata = genesets %>%
      select(-c(gene, mor, likelihood)) %>%
      distinct()
    tidy_viper_res = tdy(viper_res, "geneset", "key", "value", meta = metadata)
    return(tidy_viper_res)
  } else {
    return(viper_res)
  }
}

#' Make gene sets for VIPER
#'
#' This function convert gene sets in a table format to the format required by
#' the \code{\link[=viper]{viper::viper()}} function.
#'
#' @param genesets A dataframe of gene sets that must contain the
#' columns \code{geneset}, \code{gene}, \code{mor} and \code{likelihood}.
#'
#' @return Gene sets in the \code{\link[=viper]{viper::viper()}} format.
#'
#' @export
#' @importFrom stats setNames
make_viper_genesets = function(genesets) {
  genesets %>%
    split(.$geneset) %>%
    map(function(gs) {
      targets = setNames(gs$mor, gs$gene)
      likelihood = gs$likelihood
      list(tfmode = targets, likelihood = likelihood)
    })
}

#' Helper function
#'
#' Helper function to convert DoRothEA gene sets to standardized gene sets
#' suitable for the \code{\link[=viper]{viper::viper()}} function.
#'
#' @param genesets dorothea gene sets
#'
#' @return dorothea gene sets suitable for viper statistic
#' @export
dorothea2viper = function(genesets) {
  genesets %>%
    rename(geneset = tf, gene = target)
}

#' Helper function
#'
#' Helper function to convert PROGENy gene sets to standardized gene sets
#' suitable for the \code{\link[=viper]{viper::viper()}} function.
#'
#' @param genesets progeny gene sets
#'
#' @return progeny gene sets suitable for viper statistic
#'
#' @export
progeny2viper = function(genesets) {
  genesets %>%
    rename(geneset = pathway) %>%
    mutate(mor = sign(weight),
           likelihood = 1) %>%
    select(-weight)
}


