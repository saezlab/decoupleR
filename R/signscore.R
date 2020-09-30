##### Functions  -----------

#' SingScore Wrapper
#'
#' @param expr_matrix gene expression matrix
#' @param tiesMethod methods to use when ranking ties: min, max, average
#' @param singscore_genesets gene sets
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used for the calculation
#' @param directed
#' @param permutate Bool: Used to indicate whether we wish to obtain p-values or
#' just singscores per sample
#' @return a dataframe with pvalues for each sample or a list of singscore
#' dataframes (if permute is TRUE)
#'
#' @export
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import purrr
#' @import tidyr
#' @import singscore
run_singscore = function(expr_matrix,
                         tiesMethod,
                         singscore_genesets,
                         perm_num,
                         ncores,
                         directed,
                         permutate) {


  # Rank genes by the gene expression intensities
  rankData <- rankGenes(expr_matrix, tiesMethod = tiesMethod)

  # Go through all TFs and assign to DF   -------------
  if(directed){
    tfs_up <- names(singscore_genesets$genesets_up)
    tfs_dn <- names(singscore_genesets$genesets_dn)
    tfs_all <- union(tfs_up,tfs_dn)
    tfs_bidirectional <- intersect(tfs_up,tfs_dn)
  } else {
    tfs_all <- names(singscore_genesets)
    tfs_bidirectional <- NULL
  }



  if(permutate){
    signscore_output <- matrix(ncol = ncol(expr_matrix), nrow = 0)
    colnames(signscore_output) <- colnames(expr_matrix)
  } else{
    signscore_output <- list()
  }

  i = 1
  for(geneset in tfs_all){
    if(geneset %in% tfs_bidirectional){ ## check if gene set is bidir
      results_signscore <- simpleScore_bidir(singscore_genesets,
                                             geneset,
                                             rankData,
                                             perm_num,
                                             ncores,
                                             permutate)
    } else {
      results_signscore <- simpleScore_onedir(singscore_genesets,
                                              geneset,
                                              rankData,
                                              perm_num,
                                              ncores,
                                              directed,
                                              permutate
                                              )
    }

    # Formulate results
    if(permutate){
      signscore_output <- rbind(signscore_output, results_signscore$pvals)
      rownames(signscore_output)[i] <- geneset
    } else {
      signscore_output[[i]] <- results_signscore
      names(signscore_output)[i] <- geneset
    }
    print(paste(i, "out of", length(tfs_all), sep = " "))
    i=i+1
  }
  return(signscore_output)
}




#' Wrapper function to calculate singScore and pvalues for bidirectional sets
#' which implements the \code{\link[=singscore]{singscore::simpleScore()}} and
#' implements the \code{\link[=singscore]{singscore::getPvals()}} functions.
#'
#' @param rankData ranked expr matrix
#' @param singscore_genesets dorothea regulon gene sets with both positively
#'  and negatively regulated genes
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
simpleScore_bidir = function(singscore_genesets, geneset, rankData,
                             perm_num, ncores, permutate) {
  # Get positively and negatively regulated genes
  dorothea_up <- singscore_genesets$genesets_up[[geneset]]
  dorothea_dn <- singscore_genesets$genesets_dn[[geneset]]

  # Calculate singscores
  scoredf <- simpleScore(rankData, upSet = dorothea_up, downSet = dorothea_dn)

  if(permutate){
    permuteResult <- generateNull_bidirect(dorothea_up,
                                           dorothea_dn,
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


#' Function to calculate singScore and pvalues for onedirectional sets
#' which implements the \code{\link[=singscore]{singscore::simpleScore()}} and
#' implements the \code{\link[=singscore]{singscore::getPvals()}} functions.
#'
#' @param rankData ranked expr matrix
#' @param singscore_genesets dorothea regulon gene sets
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
simpleScore_onedir = function(singscore_genesets, geneset, rankData,
                              perm_num, ncores, directed, permutate) {


  if(directed){
    tfs_up <- names(singscore_genesets$genesets_up)
    tfs_dn <- names(singscore_genesets$genesets_dn)

    # check if the TF is in the up-regulated set
    if(geneset %in% setdiff(tfs_up,tfs_dn)){
      geneset_members <- singscore_genesets$genesets_up[[geneset]]
    } else{
      geneset_members <- singscore_genesets$genesets_dn[[geneset]]
    }

  } else{
    geneset_members <- singscore_genesets[[geneset]]
    # remove duplicates
    geneset_members <- geneset_members[!duplicated(geneset_members)]
  }

  # calculate singscore
  scoredf <- simpleScore(rankData,
                         upSet = geneset_members,
                         knownDirection = FALSE)

  if(permutate){
    permuteResult <- generateNull_one_direct(geneset_members,
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
#' of signscore will be compared (One directional - TFs with + and - signs)
#'
#' Refer to the \code{\link[=singscore]{singscore::generateNull()}} function.
#'
#' @param dorothea_up A character vector of positively reg. gene set IDs
#' @param dorothea_dn A character vector of negatively reg. gene set IDs
#' @param rankData Ranked expr_matrix
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used for the calculation
#'
#' @return permutation results to be used to perform a one-tailed test
#' @keywords internal
generateNull_bidirect = function(dorothea_up,
                                 dorothea_dn,
                                 rankData,
                                 perm_num,
                                 ncores) {
  permuteResult <-
    generateNull(
      upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
      downSet = dorothea_dn, # down-regulated gene set
      rankData = rankData,
      subSamples = 1:ncol(rankData), # which samples
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
#' of signscore will be compared (Bidirectional - TFs with + OR - signs)
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
generateNull_one_direct = function(geneset_members, rankData, perm_num, ncores){
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
#' Helper function to convert DoRothEA gene sets to two seperate sets of
#' TFs with positively and negatively regulated targets
#' suitable for the \code{\link[=singscore]{singscore::simpleScore()}} function.
#'
#' @param regulons dorothea regulon gene sets
#'
#' @return positively and negatively regulated dorothea gene sets suitable for
#' signscore
#' @keywords internal
dorothea2singscore = function(regulons) {
  # get positively and negatively regulated gene sets
  regulons_up <- regulons[regulons$mor==1,]
  regulons_dn <- regulons[regulons$mor==-1,]

  # Positively regulated genes
  genesets_up <- list()
  i = 1
  targets <- c()
  for(geneset in unique(regulons_up$geneset)){
    targets <- regulons_up[regulons_up$geneset==geneset,][["gene"]]
    genesets_up[[i]] <- targets
    names(genesets_up)[i] <- geneset
    i=i+1
    targets <- c()
  }

  # Negatively regulated genes
  genesets_dn <- list()
  i = 1
  targets <- c()
  for(geneset in unique(regulons_dn$geneset)){
    targets <- regulons_dn[regulons_dn$geneset==geneset,][["gene"]]
    genesets_dn[[i]] <- targets
    names(genesets_dn)[i] <- geneset
    i=i+1
    targets <- c()
  }

  genesets <- list(genesets_up, genesets_dn)
  names(genesets) <- c("genesets_up", "genesets_dn")

  return(genesets)
}




#' Helper function
#'
#' Helper function to convert RegNetwork gene sets to singscore format
#' suitable for the \code{\link[=singscore]{singscore::simpleScore()}} function.
#'
#' @param regulons regnetowrk regulon gene sets
#'
#' @return a named list of lists where each gene set is a list
#' @keywords internal
regnetwork2singscore = function(regulons) {

  # Undirected gene regulation
  genesets <- list()
  targets <- c()
  i = 1
  for(geneset in unique(regulons$geneset)){
    targets <- regulons[regulons$geneset==geneset,][["gene"]]
    genesets[[i]] <- targets
    names(genesets)[i] <- geneset
    i=i+1
    targets <- c()
  }

  return(genesets)
}

