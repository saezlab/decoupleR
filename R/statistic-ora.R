#' Over Representation Analysis (ORA)
#'
#' @description
#' Calculates regulatory activities using ORA.
#'
#' @details
#' ORA measures the overlap between the target feature set and a list of most
#' altered molecular features in mat. The most altered molecular features can
#' be selected from the top and or bottom of the molecular readout distribution,
#' by default it is the top 5% positive values. With these, a contingency table
#' is build and a one-tailed Fisher’s exact test is computed to determine if a
#' regulator’s set of features are over-represented in the selected features
#' from the data. The resulting score, `ora`, is the minus log10 of the
#' obtained p-value.
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
#' @param seed A single value, interpreted as an integer, or NULL for random
#'  number generation.
#' @param minsize Integer indicating the minimum number of targets per source.
#' @inheritDotParams stats::fisher.test -x -y
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `statistic`: Indicates which method is associated with which score.
#'  2. `source`: Source nodes of `network`.
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
#' run_ora(mat, network, .source='tf')
run_ora <- function(mat,
                    network,
                    .source = .data$source,
                    .target = .data$target,
                    n_up = ceiling(0.05 * nrow(mat)),
                    n_bottom = 0,
                    n_background = 20000,
                    with_ties = TRUE,
                    seed = 42,
                    minsize = 5,
                    ...) {
    # Check for NAs/Infs in mat
    check_nas_infs(mat)

    # Before to start ---------------------------------------------------------
    network <- network %>%
      rename_net({{ .source }}, {{ .target }})
    network <- filt_minsize(rownames(mat), network, minsize)
    regulons <- extract_sets(network)

    ns <- .ora_check_ns(n_up, n_bottom, n_background, network, mat)
    n_up <- ns[1]
    n_bottom <- ns[2]
    n_background <- ns[3]

    withr::with_seed(seed, {
      targets <- .ora_slice_targets(mat, n_up, n_bottom, with_ties)
    })

    # Run analysis ------------------------------------------------------------
    .ora_analysis(regulons, targets, n_background, ...)
}

# Helper functions --------------------------------------------------------
#' Wrapper to execute `run_ora()` logic one finished preprocessing of data
#'
#' @inheritParams run_ora
#' @param regulons Named list; names from `.data$source` and values
#'  from `.data$target`.
#' @param targets Named list; names from columns of `mat` and
#'  values from sliced data of `mat`.
#'
#' @inherit run_scira return
#' @keywords internal
#' @noRd
.ora_analysis <- function(regulons, targets, n_background, ...) {
    expand_grid(source = names(regulons), condition = names(targets)) %>%
        rowwise(.data$source, .data$condition) %>%
        summarise(.ora_fisher_exact_test(
            expected = regulons[[.data$source]],
            observed = targets[[.data$condition]],
            n_background = n_background,
            ...
        ),
        .groups = "drop"
        ) %>%
        select(.data$source, .data$condition,
               p_value = .data$p.value, everything()
        ) %>%
        mutate(score = -log10(.data$p_value)) %>%
        add_column(statistic = "ora", .before = 1) %>%
        select(.data$statistic, .data$source, .data$condition, .data$score, .data$p_value)
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
    false_positive <- setdiff(expected, observed) %>% length()
    false_negative <- setdiff(observed, expected) %>% length()
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
.ora_slice_targets <- function(mat, n_up, n_bottom, with_ties) {
  mat %>%
    as_tibble(rownames = "target") %>%
    tidyr::pivot_longer(
      cols = -.data$target,
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(rand=stats::rnorm(n())) %>% 
    arrange(.data$condition, .data$value, .data$rand) %>%
    group_by(.data$condition) %>%
    {
      bind_rows(
        slice_max(., .data$value, n = n_up, with_ties = with_ties),
        slice_min(., .data$value, n = n_bottom, with_ties = with_ties)
      )
    } %>%
    arrange(.data$condition) %>%
    summarise(
      targets = rlang::set_names(list(.data$target), .data$condition[1]),
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
.ora_check_ns <- function(n_up, n_bottom, n_background, network, mat) {
    if (is.null(n_background)) {
        n_background <- network %>%
            pull(.data$target) %>%
            unique() %>%
            union(rownames(mat)) %>%
            length()
    } else if (n_background < 0) {
        abort("`n` must be a non-missing positive number.")
    }

    if (n_up + n_bottom >= nrow(mat)) {
        n_up <- nrow(mat)
        n_bottom <- 0
    }

    c(n_up, n_bottom, n_background)
}
