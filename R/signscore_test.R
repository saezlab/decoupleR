
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













###### V. Test with RegNetwork

# transform RegNetwork

# Load Prerequisites
setwd("~/Repos/decoupleR")

library(tidyverse)
library(here)
library(UpSetR)
library(cowplot)
library(singscore)
library(GSEABase)

###### Dorothea
# transform dorothea data  ---------------------------
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C")) %>%
  dplyr::rename(geneset = tf, gene = target)
head(regulons)

dorothea_genesets <- dorothea2singscore(regulons, minsize = 4)
head(dorothea_genesets)




###### RegNetwork --------
# transform RegNetwork
# regnetwork_high <- read.csv("./inst/testdata/inputs/high_confidence.csv")
# regnetwork_med <- read.csv("./inst/testdata/inputs/medium_confidence.csv")
# regnetwork_low <- read.csv("./inst/testdata/inputs/low_confidence.csv")
#
# regnetwork <- rbind(regnetwork_low, regnetwork_med, regnetwork_high)
# head(regnetwork)
#
# regulons <- regnetwork %>%
#   dplyr::select(regulator_symbol, target_symbol, confidence) %>%
#   dplyr::rename(geneset = regulator_symbol, gene = target_symbol) %>%
#   dplyr::filter(!grepl("miR",gene)) %>%
#   dplyr::filter(!grepl("miR",geneset)) %>%
#   dplyr::filter(confidence %in% c("High")) %>%
#   as.tibble(regnetwork)
# head(regulons)






# Filtered data by Celina
regnetwork <- read.csv("./inst/testdata/inputs/regnetwork_filtered.csv")
regulons <- regnetwork %>%
  dplyr::rename(geneset = tf, gene = target) %>%
  dplyr::filter(!grepl("miR",gene)) %>%
  dplyr::filter(!grepl("miR",geneset)) %>%
  as.tibble(regnetwork)
head(regulons)



regnetwork_genesets <- regnetwork2singscore(regulons, 4)
head(regnetwork_genesets)





#### -------
# load data
expr_matrix <- readRDS("./inst/testdata/inputs/input-expr_matrix.rds")


#### Run signscore with Dorothea
singscore_dorothea <-
  run_singscore(expr_matrix, "min", dorothea_genesets, 1000, 6, TRUE, TRUE)
head(singscore_dorothea)


#### Run signscore with Regnetwork
singscore_regnetwork <-
  run_singscore(expr_matrix, "min", regnetwork_genesets, 1000, 6, FALSE, TRUE)
head(singscore_regnetwork)


#### Run signscore with PROGENy

progeny_data <- readRDS("./inst/testdata/inputs/input-progeny_genesets.rds")
head(progeny_data)

progeny_geneset <- progeny_data %>%
  dplyr::rename(geneset = pathway) %>%
  dplyr::mutate(mor = sign(weight),
         likelihood = 1) %>%
  dplyr::select(-weight) %>%
  dplyr::select(-likelihood)

progeny_genesets <- dorothea2singscore(progeny_geneset, 4)

singscore_progeny <-
  run_singscore(expr_matrix, "min", progeny_genesets, 1000, 6, TRUE, TRUE)
head(singscore_progeny)





# transform RegNetwork

# Load Prerequisites
setwd("~/Repos/decoupleR")

library(tidyverse)
library(here)
library(UpSetR)
library(cowplot)
library(singscore)
library(GSEABase)

###### Dorothea
# transform dorothea data  ---------------------------
data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B","C")) %>%
  dplyr::rename(geneset = tf, gene = target)
head(regulons)

dorothea_genesets <- dorothea2singscore(regulons, minsize = 4)
head(dorothea_genesets)




###### RegNetwork --------



# Filtered data by Celina
regnetwork <- read.csv("./inst/testdata/inputs/regnetwork_filtered.csv")
regulons <- regnetwork %>%
  dplyr::rename(geneset = tf, gene = target) %>%
  dplyr::filter(!grepl("miR",gene)) %>%
  dplyr::filter(!grepl("miR",geneset)) %>%
  as.tibble(regnetwork)
head(regulons)



regnetwork_genesets <- regnetwork2singscore(regulons, 4)
head(regnetwork_genesets)





#### -------
# load data
expr_matrix <- readRDS("./inst/testdata/inputs/input-expr_matrix.rds")


#### Run signscore with Dorothea
singscore_dorothea <-
  run_singscore(expr_matrix, "min", dorothea_genesets, 1000, 6, TRUE, TRUE)
head(singscore_dorothea)


#### Run signscore with Regnetwork
singscore_regnetwork <-
  run_singscore(expr_matrix, "min", regnetwork_genesets, 1000, 6, FALSE, TRUE)
head(singscore_regnetwork)


#### Run signscore with PROGENy

progeny_data <- readRDS("./inst/testdata/inputs/input-progeny_genesets.rds")
head(progeny_data)

progeny_geneset <- progeny_data %>%
  dplyr::rename(geneset = pathway) %>%
  dplyr::mutate(mor = sign(weight),
                likelihood = 1) %>%
  dplyr::select(-weight) %>%
  dplyr::select(-likelihood)

progeny_genesets <- dorothea2singscore(progeny_geneset, 4)

singscore_progeny <-
  run_singscore(expr_matrix, "min", progeny_genesets, 1000, 6, TRUE, TRUE)
head(singscore_progeny)






#-------------------------------------------------------------------------
# GenesetCollection + MultiScoe



gs1 <- GeneSet(setName = names(dorothea_genesets$genesets_up)[1],
               geneIds = dorothea_genesets$genesets_up[[1]])

gs2 <- GeneSet(setName = names(dorothea_genesets$genesets_up)[2],
               geneIds = dorothea_genesets$genesets_up[[2]])

gs3 <- GeneSet(setName = "ATF3",
               geneIds = "empty")

gs4 <- GeneSet(setName = names(dorothea_genesets$genesets_dn)[1],
               geneIds = dorothea_genesets$genesets_dn[[1]])

gs5 <- GeneSet(setName = names(dorothea_genesets$genesets_dn)[2],
               geneIds = dorothea_genesets$genesets_dn[[2]])

gs6 <- GeneSet(setName = "AHR",
               geneIds = "empty")

gsc_up <- GeneSetCollection(gs1, gs2, gs3)
gsc_dn<- GeneSetCollection(gs6, gs4, gs5)
gsc_up
gsc_dn


# Rank genes by the gene expression intensities
rankData <- rankGenes(expr_matrix, tiesMethod = "min")
head(rankData)

# remove duplicates
# rankData <- rankData [!duplicated(rankData, by = "row.id"),]



multi_scoredf <- multiScore(rankData,
                            gsc_dn,
                            subSamples = NULL,
                            centerScore = FALSE,
                            dispersionFun = mad,
                            knownDirection = FALSE)

multi_scoredf


GeneSetCollectiongsc_up[2]
geneIds(gsc_up[2])




generateNull_one_direct(gsc_up, rankData, 100, 6)

gsc_up

head(permuteResult)

# Estimate empirical p-values
pvals <- getPvals(permuteResult, scoredf, subSamples = 1:ncol(rankData))
pvals





#------------------------------
# SingScore Final Test runs


# Load Prerequisites
setwd("~/Repos/decoupleR")

library(tidyverse)
library(here)
library(UpSetR)
library(cowplot)
library(singscore)
library(GSEABase)

###### Dorothea
# transform dorothea data  ---------------------------
data(dorothea_hs, package = "dorothea")

dorothea <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A","B","C"))

genesets <- dorothea2viper(dorothea)

dorothea_genesets <- directed2singscore(genesets, minsize = 4)
head(dorothea_genesets)


# RegNetwork filtered data by Celina
regnetwork <- read.csv("./inst/testdata/inputs/regnetwork_filtered.csv")
genesets <- regnetwork2singscore(regnetwork)


regnetwork_genesets <- undirected2singscore(genesets, 4)
head(regnetwork_genesets)





#### -------
# load data
expr_matrix <- readRDS("./inst/testdata/inputs/input-expr_matrix.rds")


#### Run signscore with Dorothea
singscore_dorothea <-
  run_singscore(expr_matrix, "min", dorothea_genesets, 1000, 6, TRUE, FALSE)
head(singscore_dorothea)


#### Run signscore with Regnetwork
singscore_regnetwork <-
  run_singscore(expr_matrix, "min", regnetwork_genesets, 1000, 6, FALSE, FALSE)
head(singscore_regnetwork)


#### Run signscore with PROGENy

progeny_data <- readRDS("./inst/testdata/inputs/input-progeny_genesets.rds")
head(progeny_data)

progeny_data %>%

progeny <- progeny2viper(progeny_data)




progeny_singscore <- directed2singscore(progeny, 4)

singscore_progeny <-
  run_singscore(expr_matrix, "min", progeny_singscore, 1000, 6, TRUE, TRUE)
head(singscore_progeny)



####################

# Load Prerequisites
setwd("~/Repos/decoupleR")

library(tidyverse)
library(here)
library(UpSetR)
library(cowplot)
library(singscore)
library(GSEABase)

###### Dorothea
# transform dorothea data  ---------------------------
data(dorothea_hs, package = "dorothea")
dorothea <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A","B","C"))

regnetwork <- read.csv("./inst/testdata/inputs/regnetwork_filtered.csv")

progeny <- readRDS("./inst/testdata/inputs/input-progeny_genesets.rds")

expr_matrix <- readRDS("./inst/testdata/inputs/input-expr_matrix.rds")


# RUN -----

dorothea_singscore <- run_singscore(expr_matrix,
                                    "min",
                                    dorothea,
                                    "dorothea",
                                    4,
                                    100,
                                    6,
                                    TRUE,
                                    TRUE)

regnetwork_singscore <- run_singscore(expr_matrix,
                                      "min",
                                      regnetwork,
                                      "regnetwork",
                                      4,
                                      100,
                                      6,
                                      FALSE,
                                      TRUE)


progeny_singscore <- run_singscore(expr_matrix,
                                   "min",
                                   progeny,
                                   "progeny",
                                   4,
                                   100,
                                   6,
                                   TRUE,
                                   TRUE)

