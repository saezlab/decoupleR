library(dplyr)
library(stringr)
library(purrr)
library(decoupleR)

# Data to run -------------------------------------------------------------

mat <- readRDS(
  system.file("testdata/inputs", "input-expr_matrix.rds", package = "decoupleR")
)

dorothea_genesets <- readRDS(
  system.file("testdata/inputs", "input-dorothea_genesets.rds", package = "decoupleR")
)

statistics <- c(
  "scira",
  "pscira",
  "mean"
)

partial_decouple <- partial(
  .f = decouple,
  mat = mat,
  network = dorothea_genesets,
  .source = tf,
  .target = target
)

# decouple statistics -----------------------------------------------------

test_that("decouple same results as independent functions", {

  # n-statistics against no options.
  res_1 <- partial_decouple(
    .options = list(),
    statistics = statistics
  )

  # n-statistics against 1-option.
  res_2 <- partial_decouple(
    .options = list(.mor = "mor"),
    statistics = statistics
  )

  # n-statistics against n-options.
  res_3 <- partial_decouple(
    .options = list(
      scira = list(.mor = "mor"),
      pscira = list(.mor = "mor"),
      mean = list(.mor = "mor")
    ),
    statistics = statistics
  )

  expect_equal(res_1, res_2)
  expect_equal(res_2, res_3)

  # Compare results.
  res_4 <- res_1 %>%
    select(-run_id) %>%
    arrange(statistic, tf, condition)

  default_dorothea_statistics <- list.files(
    path = system.file("testdata/outputs/", package = "decoupleR"),
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

test_that("decouple tidy rename", {

  # 1-statistic against n-options.
  # Since we expect all conversion functions to use the same
  # rename function (i.e convert_f_defaults()), then this tests
  # for all methods the different tidy selection rename.
  walk(statistics, ~ {

    res_x <- partial_decouple(
      .options = list(
        string = list(.mor = "mor"),
        sym = list(.mor = sym("mor")),
        quo = list(.mor = quo(mor)),
        position = list(.mor = 4)
      ),
      statistics = c(.x)
    ) %>%
      select(-.data$run_id) %>%
      distinct()

    exp_x <- readRDS(
      system.file(str_glue("testdata/outputs/{.x}/output-{.x}_dorothea_default.rds"),
        package = "decoupleR"
      )
    ) %>%
      select(statistic, tf, condition, score, everything())

    expect_equal(res_x, exp_x)
  })
})
