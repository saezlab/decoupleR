library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", "viper", package = "decoupleR")

# Data to run -------------------------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()


# test for run_viper() ----------------------------------------------------

test_that("test run_viper with dorothea gene sets", {
    res_1 <- run_viper(emat, dorothea_genesets, verbose = FALSE) %>%
        select(-statistic_time)
    exp_1 <- file.path(expected_dir, "output-viper_dorothea_default.rds") %>%
        readRDS()

    res_2 <- run_viper(emat, dorothea_genesets, tf, target, mor, likelihood, verbose = FALSE) %>%
        select(-statistic_time)
    exp_2 <- file.path(expected_dir, "output-viper_dorothea_tidy-evaluation.rds") %>%
        readRDS()

    res_3 <- run_viper(emat, dorothea_genesets, verbose = FALSE, minsize = 4) %>%
        select(-statistic_time)
    exp_3 <- file.path(expected_dir, "output-viper_dorothea_minsize4.rds") %>%
        readRDS()

    expect_equal(res_1, exp_1)
    expect_equal(res_2, exp_2)
    expect_equal(res_3, exp_3)
})


# test_that("test run_viper with progeny gene sets", {
#   res1 <- run_viper(emat, progeny_genesets, gs_resource = "progeny", tidy = F)
#   exp1 <- readRDS(
#     system.file("testdata/outputs/viper/", "output-viper_progeny_default.rds",
#       package = "decoupleR"
#     )
#   )
#
#   res2 <- run_viper(emat, progeny_genesets, gs_resource = "progeny", tidy = T)
#   exp2 <- readRDS(
#     system.file("testdata/outputs/viper/", "output-viper_progeny_default_tidy.rds",
#       package = "decoupleR"
#     )
#   )
#
#   res3 <- run_viper(emat, progeny_genesets,
#     options = list(minsize = 4),
#     gs_resource = "progeny", tidy = F
#   )
#   exp3 <- readRDS(
#     system.file("testdata/outputs/viper/", "output-viper_progeny_minsize4.rds",
#       package = "decoupleR"
#     )
#   )
#
#   # NOTE: I forced to have res_x and exp_x the same rownames as for some reason
#   # the res_x objects listed p53 and not WNT as "last" pathway.
#   expect_equal(res1[rownames(exp1), ], exp1)
#   expect_equal(res2, dplyr::arrange(exp2, key, geneset))
#   expect_equal(res3[rownames(exp3), ], exp3)
# })
