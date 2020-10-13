library(decoupleR)

emat = readRDS(
  system.file("testdata/inputs", "input-expr_matrix.rds", package = "decoupleR"))
dorothea_genesets = readRDS(
  system.file("testdata/inputs","input-dorothea_genesets.rds", package = "decoupleR"))
progeny_genesets = readRDS(
  system.file("testdata/inputs","input-progeny_genesets.rds", package = "decoupleR"))
regnetwork_genesets = readRDS(
  system.file("testdata/inputs","input-regnetwork_genesets.rds", package = "decoupleR"))

library(viper)
library(dplyr)
library(singscore)


test_that("test run_singscore with progeny gene sets", {
  res1 =  run_singscore(emat, tiesMethod="min", progeny_genesets, .source="pathway",
                         .target="gene", .target_profile="weight", minsize=0, perm=100,
                         ncores=6, directed=TRUE, tidy=FALSE)
  exp1 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_progeny_default.rds",
                package = "decoupleR")
  )

  res2 =  run_singscore(emat, tiesMethod="min", progeny_genesets, .source="pathway",
                        .target="gene", .target_profile="weight", minsize=0, perm=100,
                        ncores=6, directed=TRUE, tidy=TRUE)
  exp2 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_progeny_tidy.rds",
                package = "decoupleR")
  )

  res3 =  run_singscore(emat, tiesMethod="min", progeny_genesets, .source="pathway",
                        .target="gene", .target_profile="weight", minsize=4, perm=100,
                        ncores=6, directed=TRUE, tidy=FALSE)
  exp3 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_progeny_minsize4.rds",
                package = "decoupleR")
  )

  # NOTE: I included tolerance because of variability due to low perm number
  expect_equal(res1[rownames(exp1),], exp1, tolerance=0.2)
  expect_equal(res2, dplyr::arrange(exp2, key, pathway), tolerance=0.2)
  expect_equal(res3[rownames(exp3),], exp3, tolerance=0.2)
})



test_that("test run_singscore with dorothea genesets", {
  res1 =  run_singscore(emat, tiesMethod="min", dorothea_genesets, .source="tf",
                        .target="target", .target_profile="mor", minsize=0, perm=100,
                        ncores=6, directed=TRUE, tidy=FALSE)
  exp1 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_dorothea_default.rds",
                package = "decoupleR")
  )

  res2 =  run_singscore(emat, tiesMethod="min", dorothea_genesets, .source="tf",
                        .target="target", .target_profile="mor", minsize=0, perm=100,
                        ncores=6, directed=TRUE, tidy=TRUE)
  exp2 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_dorothea_tidy.rds",
                package = "decoupleR")
  )

  res3 =  run_singscore(emat, tiesMethod="min", dorothea_genesets, .source="tf",
                        .target="target", .target_profile="mor", minsize=4, perm=100,
                        ncores=6, directed=TRUE, tidy=FALSE)
  exp3 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_dorothea_minsize4.rds",
                package = "decoupleR")
  )

  # NOTE: I included tolerance because of variability due to low perm number
  expect_equal(res1[rownames(exp1),], exp1, tolerance=0.2)
  expect_equal(res2[rownames(exp2),], exp2, tolerance=0.2)
  expect_equal(res3[rownames(exp3),], exp3, tolerance=0.2)
})


test_that("test run_singscore with regnetwork genesets", {
  res1 =  run_singscore(emat, tiesMethod="min", regnetwork_genesets, .source="tf",
                        .target="target", .target_profile=NULL, minsize=0, perm=100,
                        ncores=6, directed=FALSE, tidy=FALSE)
  exp1 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_regnetwork_default.rds",
                package = "decoupleR")
  )

  res2 =  run_singscore(emat, tiesMethod="min", regnetwork_genesets, .source="tf",
                        .target="target", .target_profile=NULL, minsize=0, perm=100,
                        ncores=6, directed=FALSE, tidy=TRUE)
  exp2 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_regnetwork_tidy.rds",
                package = "decoupleR")
  )

  res3 =  run_singscore(emat, tiesMethod="min", regnetwork_genesets, .source="tf",
                        .target="target", .target_profile=NULL, minsize=4, perm=100,
                        ncores=6, directed=FALSE, tidy=FALSE)
  exp3 = readRDS(
    system.file("testdata/outputs/singscore/", "output-singscore_regnetwork_minsize4.rds",
                package = "decoupleR")
  )

  # NOTE: I included tolerance because of variability due to low perm number
  expect_equal(res1[rownames(exp1),], exp1, tolerance=0.2)
  expect_equal(res2[rownames(exp2),], exp2, tolerance=0.2)
  expect_equal(res3[rownames(exp3),], exp3, tolerance=0.2)
})
