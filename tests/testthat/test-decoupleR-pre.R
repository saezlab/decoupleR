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
