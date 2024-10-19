library(decoupleR)

# Directories -------------------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")


# Data to run -------------------------------------------------------------

mat <- file.path(input_dir, "mat.rds") %>%
  readRDS()

net <- file.path(input_dir, "net.rds") %>%
  readRDS()

# Test for run_wmean function ---------------------------------------------

test_that("test run_zscore", {
  res_1 <- run_zscore(mat, net, minsize=0)

  expect_equal(res_1$score[1], 3.52, tolerance=0.01)
  expect_equal(sign(range(res_1$score)), c(-1, 1))
  
  net_2 <- net
  net_2$mor <- 1
  res_2 <- run_zscore(mat, net_2, minsize=0)
  
  expect_equal(res_2$score[1], 3.89, tolerance=0.01)
  expect_equal(sign(range(res_2$score)), c(1, 1))
})
