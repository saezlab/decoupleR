library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "mlm", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
  readRDS()

net <- file.path(input_dir, "net.rds") %>%
  readRDS()

# Test for run_mlm function ---------------------------------------------

test_that("test run_mlm", {
  res_1 <- run_mlm(mat, net, minsize=0)
  exp_1 <- file.path(expected_dir, "output-mlm.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1, tolerance=0.01)
})
