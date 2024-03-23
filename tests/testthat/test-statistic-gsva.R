library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "gsva", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
    readRDS()

net <- file.path(input_dir, "net.rds") %>%
    readRDS()

# Test for run_gsva function ---------------------------------------------

test_that("test run_gsva", {
    res_1 <- run_gsva(mat, net, minsize=1L, verbose = FALSE)
    exp_1 <- file.path(expected_dir, "output-gsva.rds") %>%
        readRDS()

    expect_equal(res_1, exp_1, tolerance=0.01)
})
