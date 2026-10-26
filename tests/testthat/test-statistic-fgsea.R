library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "fgsea", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
  readRDS()

net <- file.path(input_dir, "net.rds") %>%
  readRDS()

# Test for run_fgsea function ---------------------------------------------

test_that("test run_fgsea", {
  res_1 <- run_fgsea(mat, net, minsize=0, nproc=1)
  exp_1 <- file.path(expected_dir, "output-fgsea.rds") %>%
    readRDS()

  # The p-values come from 100 permutations: they differ between fgsea
  # versions, hence we check them with an absolute tolerance
  expect_equal(
    dplyr::select(res_1, -p_value),
    dplyr::select(exp_1, -p_value),
    tolerance=0.01
  )
  expect_lt(max(abs(res_1$p_value - exp_1$p_value)), 0.25)
})
