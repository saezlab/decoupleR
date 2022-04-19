library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
    readRDS()

net <- file.path(input_dir, "net.rds") %>%
    readRDS()

# Common expressions ------------------------------------------------------

# Available statistics
statistics <- c(
    'udt',
    'mdt',
    'aucell',
    'wmean',
    'wsum',
    'ulm',
    'mlm',
    'viper',
    'gsva',
    'ora',
    'fgsea'
)

# Arguments for statistics; same order as statistics vector.
args <- list(
    udt = list(),
    mdt = list(trees=1000),
    aucell = list(nproc=1, aucMaxRank=3),
    wmean = list(),
    wsum = list(),
    ulm = list(),
    viper = list(),
    gsva = list(),
    ora = list(n_up=3, n_bottom=3),
    fgsea = list()
)

partial_decouple <- purrr::partial(
    .f = decouple,
    mat = mat,
    network = net,
    .source = source,
    .target = target,
    statistics = statistics,
    minsize = 0,
    args = args
)

# decouple() --------------------------------------------------------------

test_that("decouple same results as independent functions", {

    # Choose the same defaults as in the section on generating expected results.
    res_decouple_defaults <-  suppressWarnings(partial_decouple(
        show_toy_call = FALSE,
        include_time = TRUE
    ) %>%
        dplyr::select(-.data$run_id, -statistic_time) %>%
        dplyr::filter(statistic != 'consensus') %>%
        dplyr::arrange(.data$statistic, .data$source, .data$condition))

    exp_decouple_defaults <- file.path(
        expected_dir,
        "decouple",
        "output-decouple.rds"
    ) %>%
        readRDS() %>%
        dplyr::arrange(.data$statistic, .data$source, .data$condition)

    expect_equal(res_decouple_defaults, exp_decouple_defaults, tolerance=0.1)
})
