library(decoupleR)
library(purrr)

# Base directories definition ---------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
output_dir <- file.path("inst", "testdata", "outputs")

# Specific directories creation -------------------------------------------
# Here you need to extend the vector with the name of your new statistic.

available_statistics <- c(
    "mean",
    "pscira",
    "scira",
    "viper",
    "gsva"
)

out <- available_statistics %>%
    file.path(output_dir, .) %>%
    setNames(object = ., nm = basename(.)) %>%
    as.list()

sapply(out, dir.create, showWarnings = TRUE)

# Collect individual default outputs files to decouple() test.

out_default <- stringr::str_glue(
    "output-{available_statistics}_dorothea_default.rds"
) %>%
    map2(.x = out, .y = ., file.path)

decouple_dir <- file.path(output_dir, "decouple")
dir.create(decouple_dir, showWarnings = TRUE)

# Load data to generated outputs ------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()

progeny_genesets <- file.path(input_dir, "input-progeny_genesets.rds") %>%
    readRDS()

#----- run_viper() -------------------------------------------------------------

run_viper(emat, dorothea_genesets) %>%
    saveRDS(out_default$viper)

run_viper(emat, dorothea_genesets, tf, target, mor, likelihood) %>%
    saveRDS(file.path(out$viper, "output-viper_dorothea_tidy-evaluation.rds"))

run_viper(emat, dorothea_genesets, minsize = 4) %>%
    saveRDS(file.path(out$viper, "output-viper_dorothea_minsize4.rds"))

#----- run_scira() -------------------------------------------------------------

run_scira(emat, dorothea_genesets) %>%
    saveRDS(out_default$scira)

run_scira(emat, dorothea_genesets, tf, target, mor) %>%
    saveRDS(file.path(out$scira, "output-scira_dorothea_tidy-evaluation.rds"))

run_scira(emat, dorothea_genesets, sparse = TRUE) %>%
    saveRDS(file.path(out$scira, "output-scira_dorothea_sparse-background-calculation.rds"))

#----- run_pscira() ------------------------------------------------------------

run_pscira(emat, dorothea_genesets) %>%
    saveRDS(out_default$pscira)

run_pscira(emat, dorothea_genesets, tf, target, mor) %>%
    saveRDS(file.path(out$pscira, "output-pscira_dorothea_tidy-evaluation.rds"))

run_pscira(emat, dorothea_genesets, sparse = TRUE) %>%
    saveRDS(file.path(out$pscira, "output-pscira_dorothea_sparse-background-calculation.rds"))

#----- run_mean() -------------------------------------------------------------

run_mean(emat, dorothea_genesets, .likelihood = NULL) %>%
    saveRDS(out_default$mean)

run_mean(emat, dorothea_genesets, tf, target, mor, .likelihood = NULL) %>%
    saveRDS(file.path(out$mean, "output-mean_dorothea_tidy-evaluation.rds"))

run_mean(emat, dorothea_genesets, sparse = TRUE, .likelihood = NULL) %>%
    saveRDS(file.path(out$mean, "output-mean_dorothea_sparse-background-calculation.rds"))

#---- run_gsva() ---------------------------------------------------------------

run_gsva(emat, dorothea_genesets) %>%
    saveRDS(out_default$gsva)

run_gsva(emat, dorothea_genesets, tf, target) %>%
    saveRDS(file.path(out$gsva, "output-gsva_dorothea_tidy-evaluation.rds"))

run_gsva(emat, dorothea_genesets, min.sz = 4) %>%
    saveRDS(file.path(out$gsva, "output-gsva_dorothea_minsize4.rds"))

# decouple() --------------------------------------------------------------
# This section should be kept at the end of the file
# and should be executed every time a new statistic
# is added or any entry of the default models is modified.

map_dfr(out_default, readRDS) %>%
    dplyr::arrange(.data$statistic, .data$tf, .data$condition) %>%
    saveRDS(file.path(decouple_dir, "output-decouple_dorothea_default.rds"))
