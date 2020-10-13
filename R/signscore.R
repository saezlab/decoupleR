#' SingScore Wrapper
#'
#' @param emat Gene expression matrix
#' @param tiesMethod Methods to use when ranking ties: min, max, average
#' @param dataset Gene set resource
#' @param .source column name of geneset/pathway in dataset
#' @param .target dataset column name of target mode of regulation,s sign,etc.
#' @param .target_profile dataset column name that indicates
#'  \code{\link[=directed2singscore]{directed2singscore()}}).
#' @param minsize Minimum number of genes/targets in the gene set
#' @param perm_num Number of permutations to be performed
#' @param ncores Number of cores to be used for the calculation
#' @param directed Bool: Use to indicated whether the Data set resource counts
#' directed genesets/pathways or not
#' @param tidy Bool: Return in tidy format or not
#'
#' @return A dataframe with p-values for each sample or same thing but as a
#' cartesian product for gene set groupings (i.e. tidy format)
#'
#' @export
#' @import dplyr
#' @import tibble
#' @import stringr
#' @import purrr
#' @import tidyr
#' @importFrom rlang quo_is_missing quo_is_null
#' @import singscore
#'
#'
#' @details Please note that a MultiScore method is exists within singscore
#' which could be used to produce a DF with scores from a GeneSetCollection,
#' but this is not feasible because generateNull still has to be applied to
#' each geneset. Also, all genesets must be directed for the GeneSetCollection
#' to work with all of them or it would again require a split in un/directed gs.
#' Also, in tidy format it returns confidence for each gene set, p-values are
#' not calculated seperately for each confidence level (e.g. if the same gs has
#' low and high they will have the same p-value - check if appropriate)
run_singscore = function(emat,
                         tiesMethod,
                         dataset,
                         .source,
                         .target,
                         .target_profile = NULL,
                         minsize,
                         perm_num,
                         ncores,
                         directed,
                         tidy) {

  gs_resource <- dataset %>%
    convert_to_singscore({{ .source }}, {{ .target }}, {{ .target_profile }})

  # Rank genes by the gene expression intensities
  rankData <- rankGenes(emat, tiesMethod = tiesMethod)

  if (directed) {
    geneset_up <- names(gs_resource$genesets_up)
    geneset_dn <- names(gs_resource$genesets_dn)
    genesets_all <- union(geneset_up, geneset_dn)
    up_dn_intersect <- intersect(geneset_up, geneset_dn)
  } else {
    genesets_all <- names(gs_resource)
    up_dn_intersect <- NULL
  }

  singscore_res <-  map(genesets_all, function(geneset) {
    if (geneset %in% up_dn_intersect) {
      simpleScore_bidir(gs_resource,
                           geneset,
                           rankData,
                           perm_num,
                           ncores)
    } else {
      simpleScore_undir(gs_resource,
                           geneset,
                           rankData,
                           perm_num,
                           ncores,
                           directed)
    }
  })  %>%
    setNames(.,genesets_all) %>%
    map(., function(x) {setNames(x,NULL)}) %>%
    rbind.data.frame() %>%
    `rownames<-`(colnames(emat)) %>%
    t()


  if (tidy) {
    metadata = dataset %>%
      dplyr::select(-c({{ .target }}, {{ .target_profile }})) %>%
      {
        if("likelihood" %in% colnames(.)) select(.,-likelihood) else .
      } %>%
      distinct()
    tidy_singscore_res =
      tdy(singscore_res, {{ .source }}, "key", "value", meta = metadata)
    return(tidy_singscore_res)
  } else {
    return(singscore_res)
  }
}




#' Function to calculate singScore and pvalues for bidirectional sets
#' which implements the \code{\link[=singscore]{singscore::simpleScore()}} and
#' implements the \code{\link[=singscore]{singscore::getPvals()}} functions.
#'
#' @param rankData ranked expression matrix
#' @param gs_resource gene sets with a list of
#' named lists with positively and negatively regulated genes
#' @param geneset GO/TF/Pathway gene set
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used
#'
#' @return A matrix of p-values
#' @keywords internal
simpleScore_bidir = function(gs_resource,
                             geneset,
                             rankData,
                             perm_num,
                             ncores) {
  # Get positively and negatively regulated genes
  genesets_up <- gs_resource$genesets_up[[geneset]]
  genesets_dn <- gs_resource$genesets_dn[[geneset]]

  # Calculate singscores
  scoredf <-
    simpleScore(rankData, upSet = genesets_up, downSet = genesets_dn)

  permuteResult <- generateNull_bidirect(genesets_up,
                                         genesets_dn,
                                         rankData,
                                         perm_num,
                                         ncores)
  pvals <- getPvals(permuteResult,
                    scoredf,
                    subSamples = 1:ncol(rankData))

  return(pvals)
}


#' Function to calculate singScore and p-values for undirected sets
#' which implements the \code{\link[=singscore]{singscore::simpleScore()}} and
#' implements the \code{\link[=singscore]{singscore::getPvals()}} functions.
#'

#' @param gs_resource dorothea regulon gene sets
#' @param geneset gene set resource
#' @param rankData ranked expression matrix
#' @param perm_num number of permutations to be performed
#' @param ncores number of cores to be used
#' @param directed Bool: indicates whether the set contains signed information
#' e.g. mode of regulation/sign, etc.
#'
#' @return A matrix of p-values
#' @keywords internal
simpleScore_undir <- function(gs_resource,
                              geneset,
                              rankData,
                              perm_num,
                              ncores,
                              directed) {
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
  }

  # calculate singscore
  scoredf <- simpleScore(rankData,
                         upSet = geneset_members)

  permuteResult <- generateNull_undirected(geneset_members,
                                           rankData,
                                           perm_num,
                                           ncores)
  pvals <- getPvals(permuteResult,
                    scoredf,
                    subSamples = 1:ncol(rankData))

  return(pvals)
}



#' Function to generate null hypothesis against which the single sample scores
#' of signscore will be compared (Directed genesets with + and - signs)
#'
#' Refer to the \code{\link[=singscore]{singscore::generateNull()}} function.
#'
#' @param genesets_up A character vector of positively reg. gene set IDs
#' @param genesets_dn A character vector of negatively reg. gene set IDs
#' @param rankData Ranked expression matrix
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
#' @param rankData Ranked expression matrix
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
      B = perm_num,
      ncores = ncores,
      seed = 1,
      useBPPARAM = NULL
    )
  return(permuteResult)
}


#' Helper function
#'
#' Helper function to convert a Directed/signed gene set to two seperate
#' gene sets according to the direction in which they regulate the genes
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
    group_split(mor, .keep=FALSE) %>%
    set_names(c("genesets_dn", "genesets_up")) %>%
    map(deframe)

  return(genesets)
}



#' Helper function
#'
#' Helper function to convert undirected/unsigned resources (e.g. RegNetwork)
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
    dplyr::filter(n >= minsize) %>%
    dplyr::group_by(geneset) %>%
    dplyr::filter(!duplicated(gene)) %>%
    dplyr::summarize(tmp = list(gene)) %>%
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
  genesets <- genesets %>%
    dplyr::rename(geneset = tf, gene = target) %>%
    dplyr::filter(!grepl("miR",gene)) %>%
    dplyr::filter(!grepl("miR",geneset)) %>%
    dplyr::mutate(mor = 0) %>%
    dplyr::mutate(likelihood = 0) %>%
    as_tibble()
}
