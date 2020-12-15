# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", package = "decoupleR")

# Data to run ------------------------------------------------------------------
input_tibble <- tibble(
        set_name="dorothea",
        bench_name="dbd",
        stats_list=list(
            c("mean",
              "pscira",
              "scira",
              "viper",
              "gsva"
            )),
        opts_list=list(list( # list of options for each stat method
            scira = list(),
            pscira = list(),
            mean = list(),
            viper = list(verbose = FALSE, minsize=0),
            gsva = list(verbose = FALSE, method = "gsva", mx.diff=FALSE)
        )),
        bexpr_loc = file.path(input_dir, "input-expr_matrix.rds"),
        bmeta_loc = file.path(input_dir, "input-mat_met.rds"),
        source_loc = file.path(input_dir, "input-dorothea_genesets.rds"),
        source_col="tf",
        target_col="target",
        filter_col="confidence",
        filter_crit=list(c("A"))
    )

run_benchmark(input_tibble, .perform=FALSE)

# Test run_benchmark function --------------------------------------------------
test_that("test benchmark base format function and activity column", {
    res_1 <- run_benchmark(input_tibble, .perform = FALSE)
    exp_1 <- readRDS(file.path(expected_dir,
                               "benchmark",
                               "output-benchmark_only_format.rds"))

    res_2 <- res_1@bench_res$activity[[1]] %>% select(tf, id, statistic, score)
    exp_2 <- exp_1$activity[[1]] %>% select(tf, id, statistic, score)

    # check format function
    expect_equal(res_1@bench_res %>% select(-activity),
                 exp_1 %>% select(-activity))

    # check only the first activity tibble
    expect_equal(res_2, exp_2)
})


