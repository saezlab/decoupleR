library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "gsva", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
  readRDS()

# Test for run_gsva function ---------------------------------------------

test_that("test run_gsva with dorothea gene sets", {
  res_1 <- run_gsva(emat, dorothea_genesets, options = list(verbose = FALSE)) %>%
    select(-statistic_time)
  exp_1 <- file.path(expected_dir, "output-gsva_dorothea_default.rds") %>%
    readRDS()

  res_2 <- run_gsva(emat, dorothea_genesets, tf, target, options = list(verbose = FALSE)) %>%
    select(-statistic_time)
  exp_2 <- file.path(expected_dir, "output-gsva_dorothea_tidy-evaluation.rds") %>%
    readRDS()

  res_3 <- run_gsva(emat, dorothea_genesets, options = list(min.sz = 4, verbose = FALSE)) %>%
    select(-statistic_time)
  exp_3 <- file.path(expected_dir, "output-gsva_dorothea_minsize4.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1)
  expect_equal(res_2, exp_2)
  expect_equal(res_1, res_2)
  expect_equal(res_3, exp_3)
})
