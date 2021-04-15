library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "ora", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()

# Test for run_ora function -----------------------------------------------
test_that("test run_mora with dorothea gene sets", {
    res_1 <- run_ora(emat, dorothea_genesets)
    exp_1 <- file.path(expected_dir, "output-ora_dorothea_default.rds") %>%
        readRDS()

    expect_equal(res_1, exp_1)
    expect_error(
        object = run_ora(emat, dorothea_genesets, n_background = -1),
        regexp = "`n` must be a non-missing positive number."
    )
})
