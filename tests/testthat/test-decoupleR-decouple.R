library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
expected_dir <- system.file("testdata", "outputs", package = "decoupleR")

# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
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
    aucell = list(nproc=1),
    wmean = list(),
    wsum = list(),
    ulm = list(),
    viper = list(),
    gsva = list(),
    ora = list(n_up=300, n_bottom=300),
    fgsea = list()
)

partial_decouple <- purrr::partial(
    .f = decouple,
    mat = mat,
    network = dorothea_genesets,
    .source = tf,
    .target = target,
    statistics = statistics,
    args = args
)

# decouple() --------------------------------------------------------------

test_that("decouple same results as independent functions", {

    # Choose the same defaults as in the section on generating expected results.
    res_decouple_defaults <- partial_decouple(
        show_toy_call = FALSE,
        include_time = TRUE
    ) %>%
        dplyr::select(-.data$run_id, -statistic_time) %>%
        dplyr::filter(statistic != 'consensus') %>%
        dplyr::arrange(.data$statistic, .data$source, .data$condition)

    exp_decouple_defaults <- file.path(
        expected_dir,
        "decouple",
        "output-decouple_dorothea_default.rds"
    ) %>%
        readRDS() %>%
        dplyr::arrange(.data$statistic, .data$source, .data$condition)

    expect_equal(res_decouple_defaults, exp_decouple_defaults, tolerance=0.1)
})
