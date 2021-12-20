library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "mdt", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
  readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
  readRDS()

# Test for run_mdt function ---------------------------------------------

test_that("test run_mdt with dorothea gene sets", {
  res_1 <- run_mdt(emat, dorothea_genesets, .source='tf', trees=1000)
  exp_1 <- file.path(expected_dir, "output-mdt_dorothea_default.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1, tolerance=0.1)
})
