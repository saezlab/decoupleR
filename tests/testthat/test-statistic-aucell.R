library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "aucell", package = "decoupleR")

# Data to run -------------------------------------------------------------
mat <- file.path(input_dir, "mat.rds") %>%
  readRDS()

net <- file.path(input_dir, "net.rds") %>%
  readRDS()

# Test for run_aucell function ---------------------------------------------

test_that("test run_aucell", {
  res_1 <- run_aucell(mat, net, minsize=0, nproc=1, aucMaxRank=3)
  exp_1 <- file.path(expected_dir, "output-aucell.rds") %>%
    readRDS()

  expect_equal(res_1, exp_1, tolerance=1)
})
