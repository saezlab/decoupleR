emat = readRDS(
  system.file("testdata/inputs", "input-expr_matrix.rds",
              package = "decoupleR"))
dorothea_genesets = readRDS(
  system.file("testdata/inputs","input-dorothea_genesets.rds",
              package = "decoupleR"))
progeny_genesets = readRDS(
  system.file("testdata/inputs","input-progeny_genesets.rds",
              package = "decoupleR"))

#----- run_viper() -------------------------------------------------------------

run_viper(emat, dorothea_genesets, gs_resource="dorothea", tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/viper/output-viper_dorothea_default.rds")

run_viper(emat, dorothea_genesets, gs_resource="dorothea", tidy=TRUE) %>%
  saveRDS("inst/testdata/outputs/viper/output-viper_dorothea_default_tidy.rds")

run_viper(emat, dorothea_genesets, options = list(minsize=4), gs_resource="dorothea", tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/viper/output-viper_dorothea_minsize4.rds")

run_viper(emat, progeny_genesets, gs_resource="progeny", tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/viper/output-viper_progeny_default.rds")

run_viper(emat, progeny_genesets, gs_resource="progeny", tidy=TRUE) %>%
  saveRDS("inst/testdata/outputs/viper/output-viper_progeny_default_tidy.rds")

run_viper(emat, progeny_genesets, options = list(minsize=4), gs_resource="progeny", tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/viper/output-viper_progeny_minsize4.rds")

#----- run_scira() -------------------------------------------------------------

out_scira_dir <- file.path("inst", "testdata", "outputs", "scira")

# dorothea dataset

run_scira(emat, dorothea_genesets) %>%
  saveRDS(file.path(out_scira_dir, "output-scira_dorothea_default.rds"))

run_scira(emat, dorothea_genesets, tf, target, mor) %>%
  saveRDS(file.path(out_scira_dir, "output-scira_dorothea_tidy-evaluation.rds"))

run_scira(emat, dorothea_genesets, .sparse = TRUE) %>%
  saveRDS(file.path(out_scira_dir, "output-scira_dorothea_sparse-background-calculation.rds"))

# progeny dataset

# run_scira(emat, progeny_genesets, gene, pathway, weight) %>%
#   saveRDS(file.path(out_scira_dir, "output-scira_progeny_tidy-evaluation.rds"))
#
# run_scira(emat, progeny_genesets, gene, pathway, weight, .sparse = TRUE) %>%
#   saveRDS(file.path(out_scira_dir, "output-scira_progeny_sparse-background-calculation.rds"))

#----- run_pscira() ------------------------------------------------------------

out_pscira_dir <- file.path("inst", "testdata", "outputs", "pscira")

# dorothea dataset

run_pscira(emat, dorothea_genesets) %>%
  saveRDS(file.path(out_pscira_dir, "output-pscira_dorothea_default.rds"))

run_pscira(emat, dorothea_genesets, tf, target, mor) %>%
  saveRDS(file.path(out_pscira_dir, "output-pscira_dorothea_tidy-evaluation.rds"))

run_pscira(emat, dorothea_genesets, .sparse = TRUE) %>%
  saveRDS(file.path(out_pscira_dir, "output-pscira_dorothea_sparse-background-calculation.rds"))

#----- run_mean() ------------------------------------------------------------

out_mean_dir <- file.path("inst", "testdata", "outputs", "mean")

# dorothea dataset

run_mean(emat, dorothea_genesets, .likelihood = NULL) %>%
  saveRDS(file.path(out_mean_dir, "output-mean_dorothea_default.rds"))

run_mean(emat, dorothea_genesets, tf, target, mor, .likelihood = NULL) %>%
  saveRDS(file.path(out_mean_dir, "output-mean_dorothea_tidy-evaluation.rds"))

run_mean(emat, dorothea_genesets, sparse = TRUE, .likelihood = NULL) %>%
  saveRDS(file.path(out_mean_dir, "output-mean_dorothea_sparse-background-calculation.rds"))
