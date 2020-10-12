emat = readRDS(
  system.file("testdata/inputs", "input-expr_matrix.rds",
              package = "decoupleR"))
dorothea_genesets = readRDS(
  system.file("testdata/inputs","input-dorothea_genesets.rds",
              package = "decoupleR"))
progeny_genesets = readRDS(
  system.file("testdata/inputs","input-progeny_genesets.rds",
              package = "decoupleR"))

regnetwork_genesets = readRDS(
  system.file("testdata/inputs","input-regnetwork_genesets.rds",
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


#---- run_singscore() ----------------------------------------------------------
# progeny
run_singscore(emat, tiesMethod="min", progeny_genesets, .source="pathway",
              .target="gene", .target_profile="weight", minsize=0, perm=100,
              ncores=6, directed=TRUE, tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_progeny_default.rds")

run_singscore(emat, tiesMethod="min", progeny_genesets, .source="pathway",
              .target="gene", .target_profile="weight", minsize=0, perm=100,
              ncores=6, directed=TRUE, tidy=TRUE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_progeny_tidy.rds")

run_singscore(emat, tiesMethod="min", progeny_genesets, .source="pathway",
              .target="gene", .target_profile="weight", minsize=4, perm=100,
              ncores=6, directed=TRUE, tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_progeny_minsize4.rds")

# dorothea
run_singscore(emat, tiesMethod="min", dorothea_genesets, .source="tf",
              .target="target", .target_profile="mor", minsize=0, perm=100,
              ncores=6, directed=TRUE, tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_dorothea_default.rds")

run_singscore(emat, tiesMethod="min", dorothea_genesets, .source="tf",
              .target="target", .target_profile="mor", minsize=0, perm=100,
              ncores=6, directed=TRUE, tidy=TRUE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_dorothea_tidy.rds")

run_singscore(emat, tiesMethod="min", dorothea_genesets, .source="tf",
              .target="target", .target_profile="mor", minsize=4, perm=100,
              ncores=6, directed=TRUE, tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_dorothea_minsize4.rds")

# regnetwork
run_singscore(emat, tiesMethod="min", regnetwork_genesets, .source="tf",
              .target="target", .target_profile=NULL, minsize=0, perm=100,
              ncores=6, directed=FALSE, tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_regnetwork_default.rds")

run_singscore(emat, tiesMethod="min", regnetwork_genesets, .source="tf",
              .target="target", .target_profile=NULL, minsize=0, perm=100,
              ncores=6, directed=FALSE, tidy=TRUE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_regnetwork_tidy.rds")

run_singscore(emat, tiesMethod="min", regnetwork_genesets, .source="tf",
              .target="target", .target_profile=NULL, minsize=4, perm=100,
              ncores=6, directed=FALSE, tidy=FALSE) %>%
  saveRDS("inst/testdata/outputs/singscore/output-singscore_regnetwork_minsize4.rds")


