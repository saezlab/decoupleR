# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", package = "decoupleR")

# Zenodo Inputs
bexpr_url = "https://zenodo.org/record/4322914/files/dorothea_bench_example.rds?download=1"
bmeta_url = "https://zenodo.org/record/4322914/files/dorothea_bench_meta.rds?download=1"
source_url = "https://zenodo.org/record/4322914/files/dorothea_filtered.rds?download=1"

# input tibble with Zenodo inputs
url_tibble <- tibble(
    set_name="dorothea",
    bench_name="dbd",
    stats_list=list(c("pscira", "viper", "ora")),
    opts_list=list(list(scira = list(),
                        viper = list(verbose = FALSE, minsize=0),
                        ora = list())),
    bexpr_loc = bexpr_url,
    bmeta_loc = bmeta_url,
    source_loc = source_url,
    source_col="tf",
    target_col="target",
    filter_col="confidence",
    filter_crit=list(c("A"))
    )



# Test run_benchmark url load, performance, and downsampling -------------------
test_that("test benchmark performance evaluation and downsampling", {
    res <- run_benchmark(url_tibble,
                         .url_bool = TRUE,
                         .form = TRUE,
                         .perform = TRUE,
                        )

    res_1 <- res %>%
        pluck(., "summary", "summary_table") %>%
        select(!ends_with("time"))

    exp_1 <- readRDS(file.path(expected_dir,
                               "benchmark",
                               "output-benchmark_url_perform.rds"))

    res_2 <- res@bench_res %>%
        select(set_name, bench_name, filter_crit, statistic, activity) %>%
        mutate(prc = activity %>%
                   map(~calc_curve(df=.x,
                                   downsampling=TRUE,
                                   times=10,
                                   curve="PR"))) %>%
        pull("prc")

    exp_2 <- readRDS(file.path(expected_dir,
                               "benchmark",
                               "output-benchmark_downsample.rds"))

    # check performance summary from url load
    expect_equal(res_1, exp_1)

    # check downsampling
    expect_equal(res_2, exp_2)
})
