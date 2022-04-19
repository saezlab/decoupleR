library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "viper", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
    readRDS()

net <- file.path(input_dir, "net.rds") %>%
    readRDS()

# Test for run_viper() ----------------------------------------------------

test_that("test run_viper", {
    res_1 <- run_viper(mat, net, minsize=0, verbose = FALSE)
    exp_1 <- file.path(expected_dir, "output-viper.rds") %>%
        readRDS()

    expect_equal(res_1, exp_1, tolerance=0.01)
})
