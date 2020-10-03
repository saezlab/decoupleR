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

# Test for run_pscira function ---------------------------------------------

test_that("test run_pscira with dorothea gene sets", {
  res1 <- run_pscira(emat, dorothea_genesets)
  exp1 <- readRDS(
    system.file("testdata/outputs/pscira/", "output-pscira_dorothea_default.rds",
      package = "decoupleR"
    )
  )

  res2 <- run_pscira(emat, dorothea_genesets, tf, target, mor)
  exp2 <- readRDS(
    system.file("testdata/outputs/pscira/", "output-pscira_dorothea_tidy-evaluation.rds",
      package = "decoupleR"
    )
  )

  res3 <- run_pscira(emat, dorothea_genesets, .sparse = TRUE)
  exp3 <- readRDS(
    system.file("testdata/outputs/pscira/", "output-pscira_dorothea_sparse-background-calculation.rds",
      package = "decoupleR"
    )
  )

  expect_error(run_pscira(emat, dorothea_genesets, tff), "Can't rename columns that don't exist.")
  expect_error(run_pscira(emat, dorothea_genesets, times = 1), "Parameter 'times' must be greater than or equal to 2, but 1 was passed.")
  expect_equal(res1, exp1)
  expect_equal(res2, exp2)
  expect_equal(res3, exp3)
  expect_equal(res2, res3)
})
