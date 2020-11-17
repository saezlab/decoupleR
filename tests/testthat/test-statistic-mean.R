library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "mean", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
  readRDS()

# Test for run_mean function ---------------------------------------------

test_that("test run_mean with dorothea gene sets", {
  res_1 <- run_mean(emat, dorothea_genesets, .likelihood = NULL) %>%
    select(-statistic_time)
  exp_1 <- file.path(expected_dir, "output-mean_dorothea_default.rds") %>%
    readRDS()

  res_2 <- run_mean(emat, dorothea_genesets, tf, target, mor, .likelihood = NULL) %>%
    select(-statistic_time)
  exp_2 <- file.path(expected_dir, "output-mean_dorothea_tidy-evaluation.rds") %>%
    readRDS()

  res_3 <- run_mean(emat, dorothea_genesets, sparse = TRUE, .likelihood = NULL) %>%
    select(-statistic_time)
  exp_3 <- file.path(expected_dir, "output-mean_dorothea_sparse-background-calculation.rds") %>%
    readRDS()

  expect_error(run_mean(emat, dorothea_genesets, tff, .likelihood = NULL), regexp = "Column `tff` doesn't exist.", class = "vctrs_error_subscript_oob")
  expect_error(run_mean(emat, dorothea_genesets, times = 1, .likelihood = NULL), "Parameter 'times' must be greater than or equal to 2, but 1 was passed.")
  expect_equal(res_1, exp_1)
  expect_equal(res_2, exp_2)
  expect_equal(res_3, exp_3)
  expect_equal(res_2, res_3)
})
