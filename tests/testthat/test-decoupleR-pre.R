library(decoupleR)

test_that("test intersect_regulons", {
    inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
    
    mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
    network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
    
    f_network <- intersect_regulons(mat, network, tf, target, minsize = 5)
    expect_lt(nrow(f_network), nrow(network))
})
