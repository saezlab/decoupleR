#' SingScore Wrapper
#'
#' @param expr_matrix Gene expression matrix
#' @param tiesMethod Methods to use when ranking ties: min, max, average
#' @param genesets Gene set resource
#' @param resource_name A character indicating the gene set resource. Based on
#' this argument \code{genesets} is processed futher via diverse helper
#' functions (e.g. \code{\link[=dorothea2viper]{dorothea2viper()}} and
#'  \code{\link[=directed2singscore]{directed2singscore()}}).
#' @param minsize Minimum number of genes/targets in the gene set
#' @param perm_num Number of permutations to be performed
#' @param ncores Number of cores to be used for the calculation
#' @param directed Bool: Use to indicated whether the Dataset Resource counts
#' directed genesets/pathways or not
#' @param permutate Bool: Used to indicate whether we wish to obtain p-values or
#' just singscores per sample
#'
#' @return A dataframe with p-values for each sample or a list of singscore
#' dataframes (if permute is TRUE) or a list of dataframes with signscore and
#' dispersion (if permute is FALSE)
#'
#' @export
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import purrr
#' @import tidyr
#' @import singscore
#' @import GSEABase
#'
#' @examples
#' # Run singscore with dorothea
#' dorothea_singscore <- run_singscore(expr_matrix,
#' "min",
#' dorothea,
#' "dorothea",
#' 4,
#' 100,
#' 6,
#' TRUE,
#' TRUE)
#'
#' # Run singscore with Regnetwork
#' regnetwork_singscore <- run_singscore(expr_matrix,
#'                                       "min",
#'                                       regnetwork,
#'                                       "regnetwork",
#'                                       4,
#'                                       100,
#'                                       6,
#'                                       FALSE,
#'                                       TRUE)
#'
#' # Run singscore with PROGENy
#' progeny_singscore <- run_singscore(expr_matrix,
#'                                   "min",
#'                                    progeny,
#'                                    "progeny",
#'                                    4,
#'                                    100,
#'                                    6,
#'                                    TRUE,
#'                                    TRUE)
#'
#'
#' @details Please note that a MultiScore method is exists within singscore
#' which could be used to produce a DF with scores from a GeneSetCollection,
#' but this is not feasible because generateNull still has to be applied to
#' each geneset. Also, all genesets must be directed for the GeneSetCollection
#' to work with all of them or it would again require a split in un/directed gs
run_singscore = function(expr_matrix,
                         tiesMethod,
                         genesets,
                         resource_name,
                         minsize,
                         perm_num,
                         ncores,
                         directed,
                         permutate) {


  gs_resource = switch(resource_name,
                    "dorothea" = genesets %>%
                      dorothea2viper() %>% directed2singscore(.,minsize),
                    "progeny" = genesets %>%
                      progeny2viper() %>% directed2singscore(.,minsize),
                    "regnetwork" = genesets %>%
                      regnetwork2viper() %>% undirected2singscore(.,minsize))



  # Rank genes by the gene expression intensities
  rankData <- rankGenes(expr_matrix, tiesMethod = tiesMethod)

  # Go through all Genesets/TFs/Pathways and assign to DF   -------------
  if (directed) {
    geneset_up <- names(gs_resource$genesets_up)
    geneset_dn <- names(gs_resource$genesets_dn)
    genesets_all <- union(geneset_up, geneset_dn)
    up_dn_intersect <- intersect(geneset_up, geneset_dn)
  } else {
    genesets_all <- names(gs_resource)
    up_dn_intersect <- NULL
  }


  if (permutate) {
    signscore_output <- matrix(ncol = ncol(expr_matrix), nrow = 0)
    colnames(signscore_output) <- colnames(expr_matrix)
  } else{
    signscore_output <- list()
  }

  i = 1
  for (geneset in genesets_all) {
    if (geneset %in% up_dn_intersect) {
      ## check if gene set is bidir
      results_signscore <- simpleScore_bidir(gs_resource,
                                             geneset,
                                             rankData,
                                             perm_num,
                                             ncores,
                                             permutate)
    } else {
      results_signscore <- simpleScore_undir(gs_resource,
                                              geneset,
                                              rankData,
                                              perm_num,
                                              ncores,
                                              directed,
                                              permutate)
    }

    # Formulate results
    if (permutate) {
      signscore_output <- rbind(signscore_output, results_signscore$pvals)
      rownames(signscore_output)[i] <- geneset
    } else {
      signscore_output[[i]] <- results_signscore
      names(signscore_output)[i] <- geneset
    }

    print(paste(i, "out of", length(genesets_all), sep = " "))
    i = i + 1
  }
  return(signscore_output)
}




#' Function to calculate singScore and pvalues for bidirectional sets
#' which implements the \code{\link[=singscore]{singscore::simpleScore()}} and
#' implements the \code{\link[=singscore]{singscore::getPvals()}} functions.
#'
#' @param rankData ranked expr matrix
#' @param gs_resource gene sets with a list of
#' named lists with positively and negatively regulated genes
#' @param geneset GO/TF/Pathway gene set
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used
#'
#' @return A list with 3 data frames (if permutate is TRUE):
#' scoredf - data.frame of singscores and dispersions for each sample
#' (only scoredf is returned if permutate is FALSE)
#' permuteResult - a matrix of scored randomly generated gene sets with the same
#' number of genes as the gene set
#' pvals - vector with p-values per expr matrix sample for the given gene set
#' @keywords internal
simpleScore_bidir = function(gs_resource,
                             geneset,
                             rankData,
                             perm_num,
                             ncores,
                             permutate) {
  # Get positively and negatively regulated genes
  genesets_up <- gs_resource$genesets_up[[geneset]]
  genesets_dn <- gs_resource$genesets_dn[[geneset]]

  # Calculate singscores
  scoredf <-
    simpleScore(rankData, upSet = genesets_up, downSet = genesets_dn)

  if (permutate) {
    permuteResult <- generateNull_bidirect(genesets_up,
                                           genesets_dn,
                                           rankData,
                                           perm_num,
                                           ncores)
    pvals <- getPvals(permuteResult,
                      scoredf,
                      subSamples = 1:ncol(rankData))
    score_list <- list(scoredf, permuteResult, pvals)
    names(score_list) <- c("score", "permuteResult", "pvals")

    return(score_list)
  }

  return(scoredf)
}


#' Function to calculate singScore and p-values for undirected sets
#' which implements the \code{\link[=singscore]{singscore::simpleScore()}} and
#' implements the \code{\link[=singscore]{singscore::getPvals()}} functions.
#'
#' @param rankData ranked expr matrix
#' @param gs_resource dorothea regulon gene sets
#' @param geneset
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used
#' @param directed
#' @param permutate Bool: TRUE if we want to permute and obtain permutate,
#' FALSE if we are interested only in the scores of the dataset
#'
#' @return A list with 3 data frames (if permutate is TRUE):
#' scoredf - data.frame of singscores and dispersions for each sample
#' (only scoredf is returned if permutate is FALSE)
#' permuteResult - a matrix of scored randomly generated gene sets with the same
#' number of genes as the gene set
#' pvals - vector with p-values per expr matrix sample for the given gene set
#' @keywords internal
simpleScore_undir <- function(gs_resource,
                              geneset,
                              rankData,
                              perm_num,
                              ncores,
                              directed,
                              permutate) {
  if (directed) {
    geneset_up <- names(gs_resource$genesets_up)
    geneset_dn <- names(gs_resource$genesets_dn)

    # check if the TF is in the up-regulated set
    if (geneset %in% setdiff(geneset_up, geneset_dn)) {
      geneset_members <- gs_resource$genesets_up[[geneset]]
    } else{
      geneset_members <- gs_resource$genesets_dn[[geneset]]
    }

  } else{
    geneset_members <- gs_resource[[geneset]]
    # remove duplicates (in the case of RegNetwork)
    geneset_members <- geneset_members[!duplicated(geneset_members)]
  }

  # calculate singscore
  scoredf <- simpleScore(rankData,
                         upSet = geneset_members,
                         knownDirection = FALSE)

  if (permutate) {
    permuteResult <- generateNull_undirected(geneset_members,
                                             rankData,
                                             perm_num,
                                             ncores)
    pvals <- getPvals(permuteResult,
                      scoredf,
                      subSamples = 1:ncol(rankData))
    score_list <- list(scoredf, permuteResult, pvals)
    names(score_list) <- c("score", "permuteResult", "pvals")

    return(score_list)
  }
  return(scoredf)
}



#' Function to generate null hypothesis against which the single sample scores
#' of signscore will be compared (Directed genesets with + and - signs)
#'
#' Refer to the \code{\link[=singscore]{singscore::generateNull()}} function.
#'
#' @param genesets_up A character vector of positively reg. gene set IDs
#' @param genesets_dn A character vector of negatively reg. gene set IDs
#' @param rankData Ranked expr_matrix
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used for the calculation
#'
#' @return permutation results to be used to perform a one-tailed test
#' @keywords internal
generateNull_bidirect = function(genesets_up,
                                 genesets_dn,
                                 rankData,
                                 perm_num,
                                 ncores) {
  permuteResult <-
    generateNull(
      upSet = genesets_up,
      downSet = genesets_dn,
      rankData = rankData,
      subSamples = 1:ncol(rankData),
      centerScore = TRUE,
      knownDirection = TRUE,
      B = perm_num,
      ncores = ncores,
      seed = 1,
      useBPPARAM = NULL
    )
  return(permuteResult)
}




#' Function to generate null hypothesis against which the single sample scores
#' of signscore will be compared (Undirected genesets with + OR - signs)
#'
#' Refer to the \code{\link[=singscore]{singscore::generateNull()}} function.
#'
#' @param geneset_members A character vector of gene IDs of the gene set
#' @param rankData Ranked expr_matrix
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used for the calculation
#'
#' @return permutation results to be used to perform a one-tailed test
#' @keywords internal
generateNull_undirected = function(geneset_members,
                                   rankData,
                                   perm_num,
                                   ncores) {
  permuteResult <-
    generateNull(
      upSet = geneset_members,
      downSet = NULL,
      rankData = rankData,
      subSamples = 1:ncol(rankData),
      centerScore = FALSE,
      knownDirection = FALSE,
      B = perm_num,
      ncores = ncores,
      seed = 1,
      useBPPARAM = NULL
    )
  return(permuteResult)
}


#' Helper function
#'
#' Helper function to convert a Directed gene set to two seperate gene sets of
#' according to the direction in which they regulate the genes
#' suitable for the \code{\link[=singscore]{singscore::simpleScore()}} function.
#'
#' @param genesets Directed genesets (e.g. dorothea, progeny)
#' @param minsize Minimum number of genes/targets in the geneset
#'
#' @return A list of two lists: one for positively and one list for negatively
#' regulated genesets with gene/target character vectors suitable for signscore
#' @keywords internal
directed2singscore = function(genesets, minsize) {

  genesets <- genesets %>%
    group_by(geneset, mor) %>%
    add_count(geneset) %>%
    filter(n >= minsize) %>%
    summarize(tmp = list(gene)) %>%
    ungroup() %>%
    group_split(mor, keep=FALSE) %>%
    set_names(c("genesets_dn", "genesets_up")) %>%
    map(deframe)

  return(genesets)
}



#' Helper function
#'
#' Helper function to convert undirected resources (e.g. RegNetwork)
#' to singscore format suitable for the
#' \code{\link[=singscore]{singscore::simpleScore()}} function.
#'
#' @param genesets Undirected genesets genesets
#' @param minsize Minimum number of genes regulated
#'
#' @return a named list of lists where each element is a gene set with the
#' genes/targets within that gene set
#' @keywords internal
undirected2singscore = function(genesets, minsize) {

  genesets <- genesets %>%
    add_count(geneset) %>%
    filter(n >= minsize) %>%
    group_by(geneset) %>%
    summarize(tmp = list(gene)) %>%
    deframe()

  return(genesets)
}



#' Helper function
#'
#' Helper function to convert RegNetwork gene sets to standardized gene sets
#' \code{\link[=singscore]{singscore::simpleScore()}} function.
#'
#' @param genesets RegNetwork gene sets
#'
#' @return RegNetwork gene sets suitable for singscore statistic
#'
#' @details Please note that microRNAs were filtered when working with singscore
#' @keywords internal
regnetwork2viper = function(genesets) {
  genesets <- regnetwork %>%
    dplyr::rename(geneset = tf, gene = target) %>%
    dplyr::filter(!grepl("miR",gene)) %>%
    dplyr::filter(!grepl("miR",geneset)) %>%
    as.tibble(regnetwork)
}
