library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "mdt", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
  readRDS()

net <- file.path(input_dir, "net.rds") %>%
  readRDS()

# Test for run_mdt function ---------------------------------------------

test_that("test run_mdt", {
  res_1 <- run_mdt(mat, net, minsize=0, trees=1000)
  exp_1 <- file.path(expected_dir, "output-mdt.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1, tolerance=0.1)
})
