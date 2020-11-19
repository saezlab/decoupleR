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
  res_decouple_defaults <- decouple(
    mat = mat,
    network = dorothea_genesets,
    .source = tf,
    .target = target,
    statistics = statistics,
    .options = list(
      scira = list(),
      pscira = list(),
      mean = list(.likelihood = NULL),
      viper = list(verbose = FALSE),
      gsva = list(options = list(verbose = FALSE))
    )
  ) %>%
    dplyr::select(-.data$run_id) %>%
    dplyr::arrange(.data$statistic, .data$tf, .data$condition) %>%
    select(-.data$statistic_time)

  exp_decouple_defaults <- file.path(
    expected_dir,
    "decouple",
    "output-decouple_dorothea_default.rds"
  ) %>%
    readRDS()

  expect_equal(res_decouple_defaults, exp_decouple_defaults)
})
