library(decoupleR)

test_that("multiplication works", {
    network <- tribble(
        ~tf, ~target,
        1, 1,
        1, 2,
        1, 3,
        2, 3,
        3, 4,
        4, 5,
        4, 6
    )

    # By default does not filter any.
    network %>%
        filter_regulons(.data$tf) %>%
        expect_equal(network)
})

test_that("test intersect_regulons", {
    inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
    
    mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
    network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
    
    f_network <- intersect_regulons(mat, network, tf, target, minsize = 5)
    expect_lt(nrow(f_network), nrow(network))
})
