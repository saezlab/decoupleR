library(dplyr)
library(stringr)
library(purrr)
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

statistics <- c(
  "scira",
  "pscira",
  "mean",
  "viper",
  "gsva"
)

partial_decouple <- partial(
  .f = decouple,
  mat = mat,
  network = dorothea_genesets,
  .source = tf,
  .target = target,
  statistics = statistics
)

# decouple statistics -----------------------------------------------------

test_that("decouple same results as independent functions", {

  # n-statistics against no options.
  res_1 <- partial_decouple(
    .options = list()
  )

  # Removed since no all statistics share at least one parameter.
  # # n-statistics against 1-option.
  # res_2 <- partial_decouple(
  #   .options = list(.mor = "mor")
  # )

  # n-statistics against n-options.
  res_2 <- partial_decouple(
    .options = list(
      scira = list(.mor = "mor"),
      pscira = list(.mor = "mor"),
      mean = list(.mor = "mor"),
      viper = list(.mor = "mor"),
      gsva = list(options = list(verbose = FALSE))
    )
  )

  expect_equal(res_1, res_2)

  # Compare results.
  res_4 <- res_1 %>%
    select(-run_id) %>%
    arrange(statistic, tf, condition)

  default_dorothea_statistics <- list.files(
    path = expected_dir,
    pattern = "_dorothea_default.rds",
    full.names = TRUE,
    recursive = TRUE
  ) %>%
    str_subset(
      string = .,
      pattern = paste0(statistics, collapse = "|")
    ) %>%
    map_dfr(readRDS) %>%
    select(statistic, tf, condition, score, everything()) %>%
    arrange(statistic, tf, condition)

  expect_equal(res_4, default_dorothea_statistics)
})

# Tidy selection ----------------------------------------------------------
# This must be tested in `test-utils-dataset-converters`.
