library(decoupleR)

# Data to run -------------------------------------------------------------

emat <- readRDS(
  system.file("testdata/inputs", "input-expr_matrix.rds", package = "decoupleR")
)
dorothea_genesets <- readRDS(
  system.file("testdata/inputs", "input-dorothea_genesets.rds", package = "decoupleR")
)
progeny_genesets <- readRDS(
  system.file("testdata/inputs", "input-progeny_genesets.rds", package = "decoupleR")
)

# Test for run_scira function ---------------------------------------------

test_that("test run_scira with dorothea gene sets", {
  res1 <- run_scira(emat, dorothea_genesets)
  exp1 <- readRDS(
    system.file("testdata/outputs/scira/", "output-scira_dorothea_default.rds",
      package = "decoupleR"
    )
  )

  res2 <- run_scira(emat, dorothea_genesets, tf, target, mor)
  exp2 <- readRDS(
    system.file("testdata/outputs/scira/", "output-scira_dorothea_tidy-evaluation.rds",
      package = "decoupleR"
    )
  )

  res3 <- run_scira(emat, dorothea_genesets, .sparse = TRUE)
  exp3 <- readRDS(
    system.file("testdata/outputs/scira/", "output-scira_dorothea_sparse-background-calculation.rds",
      package = "decoupleR"
    )
  )

  expect_equal(res1, exp1)
  expect_equal(res2, exp2)
  expect_equal(res3, exp3)
  expect_equal(res2, res3)
})

test_that("test run_scira with progeny gene sets", {
  res1 <- run_scira(emat, progeny_genesets, gene, pathway, weight)
  exp1 <- readRDS(
    system.file("testdata/outputs/scira/", "output-scira_progeny_tidy-evaluation.rds",
      package = "decoupleR"
    )
  )

  res2 <- run_scira(emat, progeny_genesets, gene, pathway, weight, .sparse = TRUE)
  exp2 <- readRDS(
    system.file("testdata/outputs/scira/", "output-scira_progeny_sparse-background-calculation.rds",
      package = "decoupleR"
    )
  )

  expect_error(run_scira(emat, progeny_genesets), regexp = "Column `tf` not found in `.data`", class = "dplyr:::mutate_error")
  expect_equal(res1, exp1)
  expect_equal(res2, exp2)
  expect_equal(res1, res2)
})
