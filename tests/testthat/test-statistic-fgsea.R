library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "fgsea", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
  readRDS()

# Test for run_fgsea function ---------------------------------------------

test_that("test run_fgsea with dorothea gene sets", {
  res_1 <- run_fgsea(emat, dorothea_genesets, .source='tf', nproc=1)
  exp_1 <- file.path(expected_dir, "output-fgsea_dorothea_default.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1, tolerance=0.01)
})
