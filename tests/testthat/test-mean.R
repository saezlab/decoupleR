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

# Test for run_mean function ---------------------------------------------

test_that("test run_mean with dorothea gene sets", {
  res1 <- run_mean(emat, dorothea_genesets, .likelihood = NULL)
  exp1 <- readRDS(
    system.file("testdata/outputs/mean/", "output-mean_dorothea_default.rds",
      package = "decoupleR"
    )
  )

  res2 <- run_mean(emat, dorothea_genesets, tf, target, mor, .likelihood = NULL)
  exp2 <- readRDS(
    system.file("testdata/outputs/mean/", "output-mean_dorothea_tidy-evaluation.rds",
      package = "decoupleR"
    )
  )

  res3 <- run_mean(emat, dorothea_genesets, sparse = TRUE, .likelihood = NULL)
  exp3 <- readRDS(
    system.file("testdata/outputs/mean/", "output-mean_dorothea_sparse-background-calculation.rds",
      package = "decoupleR"
    )
  )

  expect_error(run_mean(emat, dorothea_genesets, tff, .likelihood = NULL), class = "vctrs_error_subscript_oob")
  expect_error(run_mean(emat, dorothea_genesets, times = 1, .likelihood = NULL), "Parameter 'times' must be greater than or equal to 2, but 1 was passed.")
  expect_equal(res1, exp1)
  expect_equal(res2, exp2)
  expect_equal(res3, exp3)
  expect_equal(res2, res3)
})
