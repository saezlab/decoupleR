library(decoupleR)

# Data set to test ---------------------------------------------------------

dorothea_genesets <- readRDS(
    system.file("testdata/inputs", "input-dorothea_genesets.rds", package = "decoupleR")
)

# convert_f_defaults ------------------------------------------------------

test_that("convert_f_defaults (select-transmute)-like property", {
    res_1 <- convert_f_defaults(
        .data = dorothea_genesets,
        tf = tf,
        target = target
    )

    exp_1 <- dorothea_genesets %>%
        select(tf, target)

    expect_equal(res_1, exp_1)
})

test_that("convert_f_defaults swap property for single.", {

    # Normal rename
    res_1 <- convert_f_defaults(
        .data = dorothea_genesets,
        tf = target,
        target = tf
    )

    exp_1 <- dorothea_genesets %>%
        select(tf = target, target = tf)

    expect_equal(res_1, exp_1)

    # Default dplyr::rename don't allow this.
    res_1 <- convert_f_defaults(
        .data = dorothea_genesets,
        tf = target,
        .use_dots = TRUE
    )

    exp_1 <- dorothea_genesets %>%
        select(tf = target)

    expect_equal(res_1, exp_1)

    # If use_dots its false it's like normal rename.
    expect_error(convert_f_defaults(
        .data = dorothea_genesets,
        tf = target,
        .use_dots = FALSE
    ),
    regexp = "Names must be unique.",
    class = "vctrs_error_names_must_be_unique"
    )
})

test_that("convert_f_defaults add columns with defaults", {
    res_1 <- convert_f_defaults(
        .data = dorothea_genesets,
        tf = tf,
        target = target,
        mor = NULL,
        likelihood = NULL,
        .def_col_val = c(mor = 0, likelihood = 1)
    )

    exp_1 <- dorothea_genesets %>%
        select(tf, target) %>%
        mutate(mor = 0, likelihood = 1)

    expect_equal(res_1, exp_1)

    expect_error(
        convert_f_defaults(
            .data = dorothea_genesets,
            tf = tf,
            target = target,
            .def_col_val = c(mor = 0, likelihood = 1)
        ),
        regexp = "Output columns are different than expected.\nExpected: tf, target\nExtra: mor, likelihood\nRemoved: ",
        class = "different_set_columns"
    )
})

test_that("missing quos", {
    convert_to_foo <- function(x) {
        .check_quos_status({{ x }}, .dots_names = "x")
    }

    expect_error(
        object = convert_to_foo(),
        regexp = 'Quo "x" is missing, with no default.',
        class = "quo_missing_error"
    )

    expect_error(
        object = convert_to_foo(NULL),
        regexp = 'Quo "x" can not be NULL.',
        class = "quo_null_error"
    )
})
