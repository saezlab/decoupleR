
# Load Prerequisites
setwd("~/Repos/decoupleR-benchmark")


expr_matrix <- readRDS("./inst/testdata/inputs/input-expr_matrix.rds")
tf_part_test <- readRDS("./inst/testdata/inputs/input-tf_pert_test_data.rds")
dorothea_data <- readRDS("./inst/testdata/inputs/input-dorothea_genesets.rds")
rm(tf_part_test)
rm(dorothea_data)

library(tidyverse)
library(here)
library(UpSetR)
library(cowplot)
library(singscore)
library(GSEABase)

# Example data set
colnames(tgfb_expr_10_se)

# tgfb_gs_up : up regulated gene set in the tgfb gene signature
# tgfb_gs_dn : down regulated gene set in the tgfb gene signature
tgfb_gs_up
tgfb_gs_dn

length(GSEABase::geneIds(tgfb_gs_up))
length(GSEABase::geneIds(tgfb_gs_dn))

# The recommended method for dealing with ties in ranking is 'min', you can
# change by specifying 'tiesMethod' parameter for rankGenes function.
rankData <- rankGenes(tgfb_expr_10_se)
?rankGenes


# Given the ranked data and gene signature, simpleScore returns the scores and
# dispersions for each sample

# The returned data.frame consists of the scores for the up- and down-regulated
# gene-sets along with the combined score (TotalScore)
scoredf <- simpleScore(rankData, upSet = tgfb_gs_up, downSet = tgfb_gs_dn)
scoredf

?simpleScore
head(rankData[,2,drop = FALSE])


# Estimate empirical p-values for the obtained scores in individual samples
# and plot null distributions

# Hypothesis testing of the calculated scores is performed using
# a permutation test.

# The null hypothesis is that the gene-set is not enriched in the sample.

# For each sample, gene labels are randomly shuffled
# and scores computed against the gene-set.
# This is done B times to generate the null distribution.
# The generateNull() function computes these for multiple samples (n)
# simultaneously resulting in an n×B matrix with permuted scores along
# the columns for each sample.


# Permutation test
permuteResult <-
  generateNull(
    upSet = tgfb_gs_up, # A character vector of gene IDs of up-regulated gene set
    downSet = tgfb_gs_dn, # down-regulated gene set
    rankData = rankData,
    subSamples = 1:5, # which samples
    centerScore = TRUE,
    knownDirection = TRUE,
    B = 1000,
    ncores = 1,
    seed = 1,
    useBPPARAM = NULL
  )
?generateNull

head(permuteResult)


## Estimate empirical p-values

# p-values can be estimated using the getPvals() function by providing
# the null distributions calculated above.
# Unless all permutations are exhausted (mostly infeasible),
# the minimum p-value obtained is 1B
pvals <- getPvals(permuteResult, scoredf, subSamples = 1:5)



# getPval returns p-values for each individual sample.
# show the p-values for first 5 samples
pvals

##   D_Ctrl_R1   D_TGFb_R1   D_Ctrl_R2   D_TGFb_R2 Hes_Ctrl_R1
##       0.994       0.001       0.998       0.001       0.536





###### -------------------------------------------------------------------------
# II. Get a gene set from Dorothea and apply Singscore to it
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))


# Get example set from Dorothea with up and down-reg
# think of what to do as not all sets have both up and down-reg regulons of tr
dorothea_example <- regulons[regulons$tf == "AR",]
dorothea_example



dorothea_up = dorothea_example[dorothea_example$mor>0,]
dorothea_up <- dorothea_up$target

dorothea_dn = dorothea_example[dorothea_example$mor<0,]
dorothea_dn <- dorothea_dn$target

#



# Rank genes by the gene expression intensities
rankData <- rankGenes(expr_matrix, tiesMethod = "min")
head(rankData)

# remove duplicates
# rankData <- rankData [!duplicated(rankData, by = "row.id"),]



# Permutation test
permuteResult <-
  generateNull(
    upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
    downSet = dorothea_dn, # down-regulated gene set
    rankData = rankData,
    subSamples = 1:ncol(rankData), # which samples
    centerScore = TRUE,
    knownDirection = TRUE,
    B = 1000,
    ncores = 1,
    seed = 1,
    useBPPARAM = NULL
  )

head(permuteResult)


# Estimate empirical p-values
pvals <- getPvals(permuteResult, scoredf, subSamples = 1:ncol(rankData))
pvals






# III. Apply Signscore Herpes dataset with IFNb treatment

# Load Data set
expr_matrix2 <- read.table("./inst/testdata/inputs/MockvsIFN.txt",
                           row.names = 1,
                           header = TRUE)
head(expr_matrix2)


library(org.Hs.eg.db)
geneIDs1 <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys=rownames(expr_matrix2),
                                  keytype = "ENSEMBL",
                                  columns = c("SYMBOL", "ENSEMBL"))

geneIDs1 <- subset(geneIDs1, (!duplicated(geneIDs1[["SYMBOL"]])))
geneIDs1 <- subset(geneIDs1, !is.na(geneIDs1[["SYMBOL"]]))


data_me <- merge(expr_matrix2, geneIDs1, by.x = "row.names", by.y = "ENSEMBL")

rownames(data_me) <- data_me[["SYMBOL"]]
data_me$SYMBOL <- NULL
data_me$Row.names <- NULL

head(data_me)



# Rank genes by the gene expression intensities
rankData <- rankGenes(data_me, tiesMethod = "average")
head(rankData)

# Given the ranked data and gene signature, simpleScore returns the scores and
# dispersions for each sample



# remove duplicates
# rankData <- rankData [!duplicated(rankData, by = "row.id"),]


data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))





# Get example set from Dorothea with up and down-reg
# think of what to do as not all sets have both up and down-reg regulons of tr
dorothea_example <- regulons[regulons$tf == "IRF2",]
dorothea_example


dorothea_up = dorothea_example[dorothea_example$mor>0,]
dorothea_up <- dorothea_up$target

dorothea_dn = dorothea_example[dorothea_example$mor<0,]
dorothea_dn <- dorothea_dn$target


## Two directions --------------------------------------------------------------
# Given the ranked data and gene signature, simpleScore returns the scores and
# dispersions for each sample

# The returned data.frame consists of the scores for the up- and down-regulated
# gene-sets along with the combined score (TotalScore)
scoredf <- simpleScore(rankData, upSet = dorothea_up, downSet = dorothea_dn)
scoredf
?simpleScore



# Permutation test
permuteResult <-
  generateNull(
    upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
    downSet = dorothea_dn, # down-regulated gene set
    rankData = rankData,
    subSamples = 1:ncol(rankData), # which samples
    centerScore = TRUE,
    knownDirection = TRUE,
    B = 10000,
    ncores = 4,
    seed = 1,
    useBPPARAM = NULL
  )
?generateNull


head(permuteResult)



## One directions --------------------------------------------------------------

scoredf <- simpleScore(rankData, upSet = dorothea_up,
                       knownDirection = FALSE)
scoredf
?simpleScore



# Permutation test
permuteResult <-
  generateNull(
    upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
    downSet = NULL, # down-regulated gene set
    rankData = rankData,
    subSamples = 1:ncol(rankData), # which samples
    centerScore = FALSE,
    knownDirection = FALSE,
    B = 1000,
    ncores = 4,
    seed = 1,
    useBPPARAM = NULL
  )
?generateNull

# Estimate empirical p-values
pvals <- getPvals(permuteResult, scoredf, subSamples = 1:ncol(rankData))
pvals

?getPvals




# IV. Use all TFs from dorothea to apply ---------------------------------------

# transform dorothea data  ------------------
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C"))
regulons_up <- regulons[regulons$mor==1,]
regulons_dn <- regulons[regulons$mor==-1,]


# obtain tfs
tfs_all <- unique(regulons$tf)
tfs_up <- unique(regulons_up$tf)
tfs_dn <- unique(regulons_dn$tf)


# Up-regulated genes
tf_sets_up <- list()
i = 1
targets <- c()
for(tf in tfs_up){
  targets <- regulons_up[regulons_up$tf==tf,][["target"]]
  tf_sets_up[[i]] <- targets
  names(tf_sets_up)[i] <- tf
  i=i+1
  targets <- c()
}
length(tf_sets_up)


# Inhibited genes
tf_sets_dn <- list()
i = 1
targets <- c()
for(tf in tfs_dn){
  targets <- regulons_dn[regulons_dn$tf==tf,][["target"]]
  tf_sets_dn[[i]] <- targets
  names(tf_sets_dn)[i] <- tf
  i=i+1
  targets <- c()
}
length(tf_sets_dn)


# Perform signscore ---------------------------


# Load Data set
expr_matrix2 <- read.table("./inst/testdata/inputs/MockvsIFN.txt",
                           row.names = 1,
                           header = TRUE)
head(expr_matrix2)


library(org.Hs.eg.db)
geneIDs1 <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys=rownames(expr_matrix2),
                                  keytype = "ENSEMBL",
                                  columns = c("SYMBOL", "ENSEMBL"))

geneIDs1 <- subset(geneIDs1, (!duplicated(geneIDs1[["SYMBOL"]])))
geneIDs1 <- subset(geneIDs1, !is.na(geneIDs1[["SYMBOL"]]))


data_me <- merge(expr_matrix2, geneIDs1, by.x = "row.names", by.y = "ENSEMBL")

rownames(data_me) <- data_me[["SYMBOL"]]
data_me$SYMBOL <- NULL
data_me$Row.names <- NULL

head(data_me)
rm(expr_matrix2)



# Rank genes by the gene expression intensities
rankData <- rankGenes(data_me, tiesMethod = "average")
head(rankData)



# Perform permutation
tf_of_interest <- "BCL6"


if(tf_of_interest %in% intersect(tfs_up,tfs_dn)){
  print("Directions")
  dorothea_up <- tf_sets_up[[tf_of_interest]]
  dorothea_dn <- tf_sets_dn[[tf_of_interest]]

  scoredf <- simpleScore(rankData, upSet = dorothea_up, downSet = dorothea_dn)
  scoredf

  # Permutation test
  permuteResult <-
    generateNull(
      upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
      downSet = dorothea_dn, # down-regulated gene set
      rankData = rankData,
      subSamples = 1:ncol(rankData), # which samples
      centerScore = TRUE,
      knownDirection = TRUE,
      B = 1000,
      ncores = 4,
      seed = 1,
      useBPPARAM = NULL
    )

} else {
  dorothea_up <- tf_sets_up[[tf_of_interest]]
  scoredf <- simpleScore(rankData, upSet = dorothea_up,
                         knownDirection = FALSE)
  scoredf

  # Permutation test
  permuteResult <-
    generateNull(
      upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
      downSet = NULL, # down-regulated gene set
      rankData = rankData,
      subSamples = 1:ncol(rankData), # which samples
      centerScore = FALSE,
      knownDirection = FALSE,
      B = 1000,
      ncores = 4,
      seed = 1,
      useBPPARAM = NULL
    )
  print("No Directions")
}

# Estimate empirical p-values
pvals <- getPvals(permuteResult, scoredf, subSamples = 1:ncol(rankData))
pvals
names(pvals)







#
tf_of_interest <- "BCL6"

if(tf_of_interest %in% intersect(tfs_up,tfs_dn)){ ## check if TF has directions
  print("Directions")
  dorothea_up <- tf_sets_up[[tf_of_interest]]
  dorothea_dn <- tf_sets_dn[[tf_of_interest]]

  scoredf <- simpleScore(rankData, upSet = dorothea_up, downSet = dorothea_dn)
  scoredf

  # Permutation test
  permuteResult <-
    generateNull(
      upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
      downSet = dorothea_dn, # down-regulated gene set
      rankData = rankData,
      subSamples = 1:ncol(rankData), # which samples
      centerScore = TRUE,
      knownDirection = TRUE,
      B = 1000,
      ncores = 4,
      seed = 1,
      useBPPARAM = NULL
    )

} else { ## TFs
  if(tf_of_interest %in% setdiff(tfs_up,tfs_dn)){ # TFs with o
    print("Up")
    dorothea_tf <- tf_sets_up[[tf_of_interest]]

  } else{
    dorothea_tf <- tf_sets_dn[[tf_of_interest]]
    print("Down")
  }

  scoredf <- simpleScore(rankData, upSet = dorothea_tf,
                         knownDirection = FALSE)
  scoredf

  # Permutation test
  permuteResult <-
    generateNull(
      upSet = dorothea_tf, # A character vector of gene IDs of up-reg gene set
      downSet = NULL, # down-regulated gene set
      rankData = rankData,
      subSamples = 1:ncol(rankData), # which samples
      centerScore = FALSE,
      knownDirection = FALSE,
      B = 1000,
      ncores = 4,
      seed = 1,
      useBPPARAM = NULL
    )
  print("No Directions")
}

# Estimate empirical p-values
pvals <- getPvals(permuteResult,
                  scoredf,
                  subSamples = 1:ncol(rankData))
pvals

pvalue_df <- rbind(pvalue_df, pvals)
rownames(pvalue_df)[1] <- tf_of_interest
pvalue_df




# Go through all TFs and assign to DF   -------------
pvalue_df <- matrix(ncol = ncol(data_me), nrow = 0)
colnames(pvalue_df) <- colnames(data_me)
pvalue_df



#
pvalue_df <- matrix(ncol = ncol(data_me), nrow = 0)
colnames(pvalue_df) <- colnames(data_me)
i = 1
for(tf_of_interest in tfs_all){
  if(tf_of_interest %in% intersect(tfs_up,tfs_dn)){ ## check if TF has directions
    print("Directions")
    dorothea_up <- tf_sets_up[[tf_of_interest]]
    dorothea_dn <- tf_sets_dn[[tf_of_interest]]

    scoredf <- simpleScore(rankData, upSet = dorothea_up, downSet = dorothea_dn)
    scoredf

    # Permutation test
    permuteResult <-
      generateNull(
        upSet = dorothea_up, # A character vector of gene IDs of up-regulated gene set
        downSet = dorothea_dn, # down-regulated gene set
        rankData = rankData,
        subSamples = 1:ncol(rankData), # which samples
        centerScore = TRUE,
        knownDirection = TRUE,
        B = 1000,
        ncores = 4,
        seed = 1,
        useBPPARAM = NULL
      )

  } else { ## TFs
    if(tf_of_interest %in% setdiff(tfs_up,tfs_dn)){ # TFs with o
      print("Up")
      dorothea_tf <- tf_sets_up[[tf_of_interest]]

    } else{
      dorothea_tf <- tf_sets_dn[[tf_of_interest]]
      print("Down")
    }

    scoredf <- simpleScore(rankData, upSet = dorothea_tf,
                           knownDirection = FALSE)
    scoredf

    # Permutation test
    permuteResult <-
      generateNull(
        upSet = dorothea_tf, # A character vector of gene IDs of up-reg gene set
        downSet = NULL, # down-regulated gene set
        rankData = rankData,
        subSamples = 1:ncol(rankData), # which samples
        centerScore = FALSE,
        knownDirection = FALSE,
        B = 1000,
        ncores = 4,
        seed = 1,
        useBPPARAM = NULL
      )
    print("No Directions")
  }
  pvals <- getPvals(permuteResult,
                    scoredf,
                    subSamples = 1:ncol(rankData))

  pvalue_df <- rbind(pvalue_df, pvals)
  rownames(pvalue_df)[i] <- tf_of_interest
  print(pvalue_df)
  print(nrow(pvalue_df))
  i=i+1
}

head(pvalue_df)





###### V. Turn Pipeline into functions



