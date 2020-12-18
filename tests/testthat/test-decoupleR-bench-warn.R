# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", package = "decoupleR")


# Test benchmark warnings ------------------------------------------------------
test_that("test benchmark warnings", {

    # duplicate network
    net_dup <- readRDS(file.path(input_dir,
                                 "input-dorothea_genesets.rds")) %>%
        bind_rows(., .)

    # expected duplicates
    exp_dups <- "151"

    res_dup <- readRDS(file.path(expected_dir,
                                 "benchmark",
                                 "output-benchmark_with_infs.rds"))
    exp_infs <- "20 infinite values were filtered"

    # check duplicates warning
    expect_warning(filter_sets(net_dup, "tf", NULL, NULL, 0, F), exp_dups)

    # check Infs warning
    expect_warning(bench_format(res_dup), exp_infs)

})
