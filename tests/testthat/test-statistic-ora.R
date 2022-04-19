library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "ora", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
    readRDS()

net <- file.path(input_dir, "net.rds") %>%
    readRDS()

# Test for run_ora function -----------------------------------------------
test_that("test run_ora", {
    res_1 <- run_ora(mat, net, minsize=0, n_up=3, n_bottom=3)
    exp_1 <- file.path(expected_dir, "output-ora.rds") %>%
        readRDS()

    expect_equal(res_1, exp_1, tolerance=0.01)
    expect_error(
        object = run_ora(mat, net, minsize=0, n_background = -1),
        regexp = "`n` must be a non-missing positive number."
    )
})
