#' Over Representation Analysis - Fisher Exact Test
#'
#' Performs an over-representation analysis using [stats::fisher.test()].
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param n_up Integer indicating the number of top targets to slice from mat.
#' @param n_bottom Integer indicating the number of bottom targets to slice from
#'  mat.
#' @param n_background Integer indicating the background size of the sliced
#'  targets. If not specified the number of background targets is determined by
#'  the total number of unique targets in the union of `mat` and `network`.
#' @param with_ties Should ties be kept together? The default, `TRUE`,
#'  may return more rows than you request. Use `FALSE` to ignore ties,
#'   and return the first `n` rows.
#' @inheritDotParams stats::fisher.test -x -y
#'
#' @return A long format tibble of the enrichment scores for each tf
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `tf`: Source nodes of `network`.
#'  3. `condition`: Condition representing each column of `mat`.
#'  4. `score`: Regulatory activity (enrichment score).
#' @family decoupleR statistics
#' @export
#' @examples
#' inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#' mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
#' network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
#'
#' run_ora(mat, network, tf, target)
run_ora <- function(mat,
                    network,
                    .source = .data$tf,
                    .target = .data$target,
                    thr = 0.01,
                    n_background = NULL,
                    with_ties = TRUE,
                    pval_corr = 'BH',
                    ...) {
    # Before to start ---------------------------------------------------------
    regulons <- network %>%
        convert_to_ora({{ .source }}, {{ .target }})

    ns <- .ora_check_ns(thr, n_background, network, mat)
    n_up <- ns[1]
    n_background <- ns[2]

    targets <- .ora_slice_targets(mat, n_up, with_ties)

    # Run analysis ------------------------------------------------------------
    .ora_analysis(regulons, targets, n_background, pval_corr, ...)
}

# Helper functions --------------------------------------------------------
#' Wrapper to execute `run_ora()` logic one finished preprocessing of data
#'
#' @inheritParams run_ora
#' @param regulons Named list; names from `.data$tf` and values
#'  from `.data$target`.
#' @param targets Named list; names from columns of `mat` and
#'  values from sliced data of `mat`.
#'
#' @inherit run_scira return
#' @keywords internal
#' @noRd
.ora_analysis <- function(regulons, targets, n_background, pval_corr, ...) {
    result <- expand_grid(tf = names(regulons), condition = names(targets)) %>%
        rowwise(.data$tf, .data$condition) %>%
        summarise(.ora_fisher_exact_test(
            expected = regulons[[.data$tf]],
            observed = targets[[.data$condition]],
            n_background = n_background,
            ...
        ),
        .groups = "drop"
        ) %>%
        select(.data$tf, .data$condition,
            score = .data$p.value, everything()
        ) %>%
        add_column(statistic = "ora", .before = 1)
  if (is.null(pval_corr)){
    result[['score']] <- -log10(result[['score']])
  } else{
    result[['score']] <- -log10(p.adjust(result[['score']],method = pval_corr))
  }
  result
}

#' Fisher Exact Test
#'
#' @inheritParams run_ora
#' @inheritParams .ora_contigency_table
#'
#' @return Single row summary "glance" of a object of class `htest`.
#' @keywords internal
#' @noRd
.ora_fisher_exact_test <- function(expected, observed, n_background, ...) {
    exec(
        .fn = stats::fisher.test,
        x = .ora_contingency_table(expected, observed, n_background),
        y = NULL,
        alternative='greater',
        !!!list(...)
    ) %>%
        broom::glance()
}

#' Create contingency table
#'
#' @inheritParams run_ora
#' @param expected Vector with expected targets
#' @param observed Vector with observed targets
#'
#' @return 2 x 2 matrix
#' @keywords internal
#' @noRd
.ora_contingency_table <- function(expected, observed, n_background) {
    true_positive <- intersect(observed, expected) %>% length()
    false_positive <- setdiff(observed, expected) %>% length()
    false_negative <- setdiff(expected, observed) %>% length()
    true_negative <- (n_background -
        true_positive - false_positive - false_negative)

    c(true_positive, false_positive, false_negative, true_negative) %>%
        matrix(nrow = 2, ncol = 2, byrow = FALSE)
}

#' Slice targets per condition
#'
#' @inheritParams run_ora
#' @return Named list with sliced targets per condition.
#'
#' @keywords internal
#' @noRd
.ora_slice_targets <- function(mat, n_up, with_ties) {
    mat %>%
        as_tibble(rownames = "target") %>%
        pivot_longer(
            cols = -.data$target,
            names_to = "condition",
            values_to = "value"
        ) %>%
        group_by(.data$condition) %>%
        slice_max(., abs(.data$value), n = n_up, with_ties = with_ties) %>%
        summarise(
            targets = set_names(list(.data$target), .data$condition[1]),
            .groups = "drop"
        ) %>%
        pull(.data$targets)
}

#' Check values of variables with n_prefix
#'
#' Set convenient default values for the ns so that downstream
#' functions work fine.
#'
#' @inheritParams run_ora
#'
#' @return ns modified if necessary.
#'
#' @keywords internal
#' @noRd
.ora_check_ns <- function(threshold, n_background, network, mat) {
    if (is.null(n_background)) {
        n_background <- network %>%
            pull(.data$target) %>%
            unique() %>%
            union(rownames(mat)) %>%
            length()
    } else if (n_background < 0) {
        abort("`n` must be a non-missing positive number.")
    }
    # Select % of total genes in mat
    n_up <- ceiling(threshold * nrow(mat))

    c(n_up, n_background)
}
