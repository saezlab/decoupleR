# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Zenodo Inputs
bexpr_url = "https://zenodo.org/record/4322914/files/dorothea_bench_example.rds?download=1"
bmeta_url = "https://zenodo.org/record/4322914/files/dorothea_bench_meta.rds?download=1"
source_url = "https://zenodo.org/record/4322914/files/dorothea_filtered.rds?download=1"

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

ulr_tibble <- input_tibble %>%
    mutate(filter_crit=list(c("A")),
        bexpr_loc = bexpr_url,
           bmeta_loc = bmeta_url,
           source_loc = source_url,
           stats_list=list(c("pscira","viper")),
           opts_list=list(list(scira=list(),
                               viper = list(verbose = FALSE, minsize=0))))

# Test run_benchmark function --------------------------------------------------
test_that("test run_benchmark stages (format, call to decouple, roc)", {
    res_1 <- run_benchmark(input_tibble, .perform = FALSE) %>%
        pluck(., "bench_res")
    exp_1 <- readRDS(file.path(expected_dir,
                               "benchmark",
                               "output-benchmark_only_format.rds"))

    res_2 <- res_1$activity[[1]] %>% select(tf, id, statistic, score)
    exp_2 <- exp_1$activity[[1]] %>% select(tf, id, statistic, score)

    res_3 <- run_benchmark(ulr_tibble, .url_bool = TRUE) %>%
        pluck(., "summary", "summary_table") %>%
        select(-c(statistic_time, regulon_time))
    exp_3 <- readRDS(file.path(expected_dir,
                               "benchmark",
                               "output-benchmark_url_perform.rds"))


    # check format function
    expect_equal(res_1 %>% select(-activity),
                 exp_1 %>% select(-activity))

    # check the format of only the first activity tibble
    expect_equal(res_2, exp_2)

    # check performance summary from url load
    expect_equal(res_3, exp_3)
})


