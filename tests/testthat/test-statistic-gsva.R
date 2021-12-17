library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "gsva", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()

# Test for run_gsva function ---------------------------------------------

test_that("test run_gsva with dorothea gene sets", {
    res_1 <- run_gsva(emat, dorothea_genesets, .source='tf', verbose = FALSE)
    exp_1 <- file.path(expected_dir, "output-gsva_dorothea_default.rds") %>%
        readRDS()

    expect_equal(res_1, exp_1, tolerance=0.01)
})
