library(decoupleR)

emat = readRDS(
  system.file("testdata/inputs", "input-expr_matrix.rds", package = "decoupleR"))
dorothea_genesets = readRDS(
  system.file("testdata/inputs","input-dorothea_genesets.rds", package = "decoupleR"))
progeny_genesets = readRDS(
  system.file("testdata/inputs","input-progeny_genesets.rds", package = "decoupleR"))


test_that("test run_viper with dorothea gene sets", {
  res1 = run_viper(emat, dorothea_genesets, gs_resource="dorothea", tidy=F)
  exp1 = readRDS(
    system.file("testdata/outputs/viper/", "output-viper_dorothea_default.rds",
                package = "decoupleR")
  )
  res2 = run_viper(emat, dorothea_genesets, gs_resource="dorothea", tidy=T)
  exp2 = readRDS(
    system.file("testdata/outputs/viper/", "output-viper_dorothea_default_tidy.rds",
                package = "decoupleR")
  )

  res3 = run_viper(emat, dorothea_genesets, options = list(minsize=4),
                   gs_resource="dorothea", tidy=F)
  exp3 = readRDS(
    system.file("testdata/outputs/viper/", "output-viper_dorothea_minsize4.rds",
                package = "decoupleR")
  )

  expect_equal(res1, exp1)
  expect_equal(res2, exp2)
  expect_equal(res3, exp3)
})

#
# test_that("test run_viper with progeny gene sets", {
#   res1 = run_viper(emat, progeny_genesets, gs_resource="progeny", tidy=F)
#   exp1 = readRDS(
#     system.file("testdata/outputs/viper/", "output-viper_progeny_default.rds",
#                 package = "decoupleR")
#   )
#   print(res1)
#   print(exp1)
#
#   res2 = run_viper(emat, progeny_genesets, gs_resource="progeny", tidy=T)
#   exp2 = readRDS(
#     system.file("testdata/outputs/viper/", "output-viper_progeny_default_tidy.rds",
#                 package = "decoupleR")
#   )
#
#   res3 = run_viper(emat, progeny_genesets, options = list(minsize=4),
#                    gs_resource="progeny", tidy=F)
#   exp3 = readRDS(
#     system.file("testdata/outputs/viper/", "output-viper_progeny_minsize4.rds",
#                 package = "decoupleR")
#   )
#
#   expect_equal(res1, exp1)
#   expect_equal(res2, exp2)
#   expect_equal(res3, exp3)
# })
