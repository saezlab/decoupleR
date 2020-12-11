library(decoupleR)

network <- tibble::tribble(
    ~tf, ~target, ~mor,
    1, 1, 1,
    2, 1, 0,
    2, 2, 0
)

test_that("test get_profile_of", {
    partial_get_profile_of <- purrr::partial(
        .f = get_profile_of,
        data = network,
        sources = list(tf = c(1, 2), target = c(1, 2))
    )

    expected_network <- tibble::tribble(
        ~tf, ~target, ~mor,
        1, 1, 1,
        1, 2, NA,
        2, 1, 0,
        2, 2, 0
    )

    expect_equal(partial_get_profile_of(), expected_network)

    expect_equal(
        partial_get_profile_of(values_fill = 0),
        replace_na(expected_network, list(mor = 0))
    )

    expect_equal(
        partial_get_profile_of(values_fill = list(mor = 0)),
        replace_na(expected_network, list(mor = 0))
    )
})

test_that("test pivot_wider_profile", {
    partial_pivot <- purrr::partial(
        .f = pivot_wider_profile,
        data = network,
        id_cols = tf,
        names_from = target,
        values_from = mor
    )

    df_mat <- partial_pivot(to_matrix = FALSE, to_sparse = FALSE)
    mat_mat <- partial_pivot(to_matrix = TRUE)
    sparse_mat <- partial_pivot(to_sparse = TRUE)

    expect_true(is.data.frame(df_mat))
    expect_true(is.matrix(mat_mat))
    expect_true(is(sparse_mat, "dsCMatrix") || is(sparse_mat, "dtCMatrix"))
})
