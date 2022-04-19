library(decoupleR)

test_that("test intersect_regulons", {
    inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
    
    mat <- readRDS(file.path(inputs_dir, "mat.rds"))
    net <- readRDS(file.path(inputs_dir, "net.rds"))
    
    f_net <- intersect_regulons(mat, net, source, target, minsize = 4)
    expect_lt(nrow(f_net), nrow(net))
})
