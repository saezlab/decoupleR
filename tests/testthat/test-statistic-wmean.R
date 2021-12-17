library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "wmean", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()

# Test for run_wmean function ---------------------------------------------

test_that("test run_wmean with dorothea gene sets", {
    res_1 <- run_wmean(emat, dorothea_genesets, .source='tf')
    exp_1 <- file.path(expected_dir, "output-wmean_dorothea_default.rds") %>%
        readRDS()

    expect_error(
        run_wmean(emat, dorothea_genesets, .source='tf', times = 1),
        "Parameter 'times' must be greater than or equal to 2, but 1 was passed."
    )
    expect_equal(res_1, exp_1, tolerance=0.01)
})
