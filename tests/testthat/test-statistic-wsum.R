library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "wsum", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
    readRDS()

net <- file.path(input_dir, "net.rds") %>%
    readRDS()

# Test for run_wsum function ---------------------------------------------

test_that("test run_wsum", {
    res_1 <- run_wsum(mat, net, minsize=0)
    exp_1 <- file.path(expected_dir, "output-wsum.rds") %>%
        readRDS()

    expect_error(
        run_wsum(mat, net, minsize=0, times = 1),
        "Parameter 'times' must be greater than or equal to 2, but 1 was passed."
    )
    expect_equal(res_1, exp_1, tolerance=0.01)
})
