# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", package = "decoupleR")

# Data to run ------------------------------------------------------------------
mat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()

input_tibble <- file.path(input_dir, "input-design_tibble.rds") %>%
    readRDS()

# Test run_benchmark function --------------------------------------------------
test_that("test benchmark base format function and activity column", {
    res_1 <- run_benchmark(input_tibble, .perform = F)
    exp_1 <- readRDS(file.path(expected_dir,
                               "benchmark",
                               "output-benchmark_only_format.rds"
    ))

    res_2 <- res_1@bench_res$activity[[1]] %>% select(tf, id, statistic, score)
    exp_2 <- exp_1$activity[[1]] %>% select(tf, id, statistic, score)

    # check format function
    expect_equal(res_1@bench_res %>% select(-activity),
                 exp_1 %>% select(-activity))

    # check only the first activity tibble
    expect_equal(res_2, exp_2)
})


