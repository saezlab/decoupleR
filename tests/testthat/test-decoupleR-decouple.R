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
    "scira",
    "pscira",
    "mean",
    "viper",
    "gsva"
)

# Arguments for statistics; same order as statistics vector.
args <- list(
    scira = list(),
    pscira = list(),
    mean = list(.likelihood = NULL),
    viper = list(verbose = FALSE),
    gsva = list(verbose = FALSE)
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

    # Available statistics
    statistics <- c(
        "scira",
        "pscira",
        "mean",
        "viper",
        "gsva"
    )

    # Choose the same defaults as in the section on generating expected results.
    res_decouple_defaults <- partial_decouple(
        show_toy_call = FALSE,
        include_time = TRUE
    ) %>%
        dplyr::select(-.data$run_id, -statistic_time) %>%
        dplyr::arrange(.data$statistic, .data$tf, .data$condition)

    exp_decouple_defaults <- file.path(
        expected_dir,
        "decouple",
        "output-decouple_dorothea_default.rds"
    ) %>%
        readRDS()

    expect_equal(res_decouple_defaults, exp_decouple_defaults)
})

test_that("see expected toy call", {

    # Choose the same defaults as in the section on generating expected results
    expect_snapshot(
        x = partial_decouple(show_toy_call = TRUE, include_time = FALSE)
    )
})
