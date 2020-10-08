setwd("./inst/testdata/inputs")
emat <- readRDS("input-expr_matrix.rds")
dorothea_data <- readRDS("input-dorothea_genesets.rds")
regNetwork <- read.csv("regnetwork_filtered.csv")
progeny_genesets <- readRDS("~/Repos/decoupleR/inst/testdata/inputs/input-progeny_genesets.rds")

library(viper)
library(tidyverse)
library(singscore)
library(rlang)



singscore_output <- run_singscore(emat,
                                  "min",
                                  progeny_genesets,
                                  "pathway",
                                  "gene",
                                  "weight",
                                  4,
                                  100,
                                  6,
                                  TRUE,
                                  TRUE)






singscore_output_tidy <- run_singscore(emat,
                                       "min",
                                       regNetwork,
                                       "tf",
                                       "target",
                                       NULL,
                                       4,
                                       100,
                                       6,
                                       FALSE,
                                       TRUE)


viper_output <- run_viper(emat = emat,
                          genesets = dorothea_data,
                          gs_resource = "dorothea",
                          tidy = TRUE)


viper_output2 <- run_viper(emat = emat,
                          genesets = dorothea_data,
                          gs_resource = "dorothea",
                          tidy = FALSE)
