library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
output_dir <- file.path("inst", "testdata", "outputs")

out <- list(
  mean = file.path(output_dir, "mean"),
  pscira = file.path(output_dir, "pscira"),
  scira = file.path(output_dir, "scira"),
  viper = file.path(output_dir, "viper"),
  gsva = file.path(output_dir, "gsva")
)

sapply(out, dir.create, showWarnings = TRUE)

# Load data to generated outputs ------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
  readRDS()

progeny_genesets <- file.path(input_dir, "input-progeny_genesets.rds") %>%
  readRDS()

#----- run_viper() -------------------------------------------------------------

run_viper(emat, dorothea_genesets) %>%
  saveRDS(file.path(out$viper, "output-viper_dorothea_default.rds"))

run_viper(emat, dorothea_genesets, tf, target, mor, likelihood) %>%
  saveRDS(file.path(out$viper, "output-viper_dorothea_tidy-evaluation.rds"))

run_viper(emat, dorothea_genesets, options = list(minsize = 4)
) %>%
  saveRDS(file.path(out$viper, "output-viper_dorothea_minsize4.rds"))

#----- run_scira() -------------------------------------------------------------

run_scira(emat, dorothea_genesets) %>%
  saveRDS(file.path(out$scira, "output-scira_dorothea_default.rds"))

run_scira(emat, dorothea_genesets, tf, target, mor) %>%
  saveRDS(file.path(out$scira, "output-scira_dorothea_tidy-evaluation.rds"))

run_scira(emat, dorothea_genesets, .sparse = TRUE) %>%
  saveRDS(file.path(out$scira, "output-scira_dorothea_sparse-background-calculation.rds"))

#----- run_pscira() ------------------------------------------------------------

run_pscira(emat, dorothea_genesets) %>%
  saveRDS(file.path(out$pscira, "output-pscira_dorothea_default.rds"))

run_pscira(emat, dorothea_genesets, tf, target, mor) %>%
  saveRDS(file.path(out$pscira, "output-pscira_dorothea_tidy-evaluation.rds"))

run_pscira(emat, dorothea_genesets, .sparse = TRUE) %>%
  saveRDS(file.path(out$pscira, "output-pscira_dorothea_sparse-background-calculation.rds"))

#----- run_mean() -------------------------------------------------------------

run_mean(emat, dorothea_genesets, .likelihood = NULL) %>%
  saveRDS(file.path(out$mean, "output-mean_dorothea_default.rds"))

run_mean(emat, dorothea_genesets, tf, target, mor, .likelihood = NULL) %>%
  saveRDS(file.path(out$mean, "output-mean_dorothea_tidy-evaluation.rds"))

run_mean(emat, dorothea_genesets, sparse = TRUE, .likelihood = NULL) %>%
  saveRDS(file.path(out$mean, "output-mean_dorothea_sparse-background-calculation.rds"))

#---- run_gsva() ---------------------------------------------------------------

run_gsva(emat, dorothea_genesets) %>%
  saveRDS(file.path(out$gsva, "output-gsva_dorothea_default.rds"))

run_gsva(emat, dorothea_genesets, tf, target) %>%
  saveRDS(file.path(out$gsva, "output-gsva_dorothea_tidy-evaluation.rds"))

run_gsva(emat, dorothea_genesets, options = list(min.sz = 4)) %>%
  saveRDS(file.path(out$gsva, "output-gsva_dorothea_minsize4.rds"))
