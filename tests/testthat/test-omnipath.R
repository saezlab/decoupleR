library(decoupleR)

test_that("test show_resources", {
  lst <- show_resources()
  testthat::expect_true(length(lst) > 0)
})

test_that("test get_resource", {
  df <- get_resource('TFcensus')
  testthat::expect_true(nrow(df) > 0)
})

test_that("test get_progeny human", {
  df <- get_progeny(organism = 'human')
  testthat::expect_true(nrow(df) > 0)
})

test_that("test get_progeny mouse", {
  df <- get_progeny(organism = 'mouse')
  testthat::expect_true(nrow(df) > 0)
})

test_that("test get_dorothea human", {
  df <- get_dorothea(organism = 'human')
  testthat::expect_true(nrow(df) > 0)
})

test_that("test get_dorothea mouse", {
  df <- get_dorothea(organism = 'mouse')
  testthat::expect_true(nrow(df) > 0)
})
