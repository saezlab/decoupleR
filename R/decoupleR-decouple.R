#' Evaluate multiple statistics with same input data
#'
#' Calculate the source activity per sample out of a gene expression matrix by
#' coupling a regulatory network with a variety of statistics.
#'
#' @inheritParams .decoupler_mat_format
#' @inheritParams .decoupler_network_format
#' @param statistics Statistical methods to be run sequentially. If none are 
#' provided, only top performer methods are run (mlm, ulm and wsum).
#' @param args A list of argument-lists the same length as `statistics`
#'  (or length 1). The default argument, list(NULL), will be recycled to the
#'  same length as `statistics`, and will call each function with no arguments
#'   (apart from `mat`, `network`, `.source` and, `.target`).
#' @param consensus_score Boolean whether to run a consensus score between
#' methods.
#' @param consensus_stats List of estimate names to use for the calculation 
#' of the consensus score. This is used to filter out extra estimations 
#' from some methods, for example wsum returns wsum, corr_wsum and norm_wsum. If
#'  none are provided, and also no statstics where provided, only top performer 
#'  methods are used (mlm, ulm and norm_wsum). Else, it will use all available 
#'  estimates after running all methods in the statistics argument. 
#' @param include_time Should the time per statistic evaluated be informed?
#' @param minsize Integer indicating the minimum number of targets per source.
#' @param show_toy_call The call of each statistic must be informed?
#'
#' @return A long format tibble of the enrichment scores for each source
#'  across the samples. Resulting tibble contains the following columns:
#'  1. `run_id`: Indicates the order in which the methods have been executed.
#'  2. `statistic`: Indicates which method is associated with which score.
#'  3. `source`: Source nodes of `network`.
#'  4. `condition`: Condition representing each column of `mat`.
#'  5. `score`: Regulatory activity (enrichment score).
#'  6. `statistic_time`: If requested, internal execution time indicator.
#'  7. `p_value`: p-value (if available) of the obtained score.
#' @export
#' @import purrr
#' @family decoupleR statistics
#' @examples
#' if (FALSE) {
#'     inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
#'
#'     mat <- readRDS(file.path(inputs_dir, "mat.rds"))
#'     net <- readRDS(file.path(inputs_dir, "net.rds"))
#'
#'     decouple(
#'         mat = mat,
#'         network = net,
#'         .source = "source",
#'         .target = "target",
#'         statistics = c("gsva", "wmean", "wsum", "ulm", "aucell"),
#'         args = list(
#'             gsva = list(verbose = FALSE),
#'             wmean = list(.mor = "mor", .likelihood = "likelihood"),
#'             wsum = list(.mor = "mor"),
#'             ulm = list(.mor = "mor")
#'         ),
#'         minsize = 0
#'     )
#' }
decouple <- function(mat,
                     network,
                     .source = source,
                     .target = target,
                     statistics = NULL,
                     args = list(NULL),
                     consensus_score = TRUE,
                     consensus_stats = NULL,
                     include_time = FALSE,
                     show_toy_call = FALSE,
                     minsize = 5) {

    # NSE vs. R CMD check workaround
    condition <- run_id <- score <- source <- statistic <- target <- NULL


    # If NULL use top performer methods.
    if (is.null(statistics)){
        statistics <- c('mlm','ulm','wsum')
        if (is.null(consensus_stats)) {
            consensus_stats <- c('mlm','ulm','norm_wsum')
        }
    } else if (length(statistics) == 1) {
        if (tolower(statistics)=='all') {
            statistics <- c('udt','mdt','aucell','wmean','wsum','ulm',
                            'mlm','viper','gsva','ora','fgsea')
        }
    }

    # Match statistic names with arguments
    for (stat in setdiff(statistics, names(args))) {
        args[[stat]] = list()
    }
    args <- args[names(args) %in% statistics]
    statistics <- statistics[match(names(args),statistics)]

    # Overwrite minsize
    for (name in names(args)) {
        args[[name]][['minsize']] <- minsize
    }

    # Match statistics to couple ----------------------------------------------
    statistics <- .select_statistics(statistics)

    # Evaluate statistics -----------------------------------------------------

    mat_symbol <- .label_expr({{ mat }})
    network_symbol <- .label_expr({{ network }})

    # For the moment this will only ensure that the parameters passed
    # to decoupleR are the same when invoking the functions.
    df <- map2_dfr(
        .x = statistics,
        .y = args,
        .f = .invoke_statistic,
        mat = mat,
        network = network,
        .source = {{ .source }},
        .target = {{ .target }},
        mat_symbol = {{ mat_symbol }},
        network_symbol = {{ network_symbol }},
        include_time = include_time,
        minsize = minsize,
        show_toy_call = show_toy_call,
        .id = "run_id"
    ) %>%
        select(
            run_id,
            statistic,
            source,
            condition,
            score,
            everything()
        ) %>%
        mutate(run_id = as.numeric(run_id))
    if (consensus_score){
        if (!is.null(consensus_stats)) {
            consensus <- df %>% 
                dplyr::filter(statistic %in% consensus_stats) %>%
                decoupleR::run_consensus(., include_time=include_time)
        } else {
            consensus <- decoupleR::run_consensus(df, include_time=include_time)
        }
            
        df <- dplyr::bind_rows(df, consensus)
    }
    df
}

# Helpers -----------------------------------------------------------------
#' Choose statistics to run
#'
#' It allows the user to select multiple statistics to run,
#' no matter if they are repeated or not.
#'
#' @details
#' From the user perspective, this could be useful since any traceback
#' would look something like decoupleR::run_{statistic}().
#'
#' @inheritParams decouple
#'
#' @return list of expressions of statistics to run.
#' @keywords internal
#' @noRd
.select_statistics <- function(statistics) {
    available_statistics <- list(
        aucell = expr(run_aucell),
        udt = expr(run_udt),
        mdt = expr(run_mdt),
        wmean = expr(run_wmean),
        ulm = expr(run_ulm),
        mlm = expr(run_mlm),
        wsum = expr(run_wsum),
        viper = expr(run_viper),
        gsva = expr(run_gsva),
        ora = expr(run_ora),
        fgsea = expr(run_fgsea)
    )

    statistics %>%
        match.arg(names(available_statistics), several.ok = TRUE) %>%
        available_statistics[.] %>%
        unname()
}

#' Construct an expression to evaluate a decoupleR statistic.
#'
#' @details
#' `.invoke_statistic()` was designed because [purrr::invoke_map_dfr()] is
#' retired. The alternative proposed by the developers by purrr is to use
#' [rlang::exec()] in combination with [purrr::map2()], however, the function
#' is not a quoting function, so the parameters that require the
#' `curly-curly` (`{{}}`) operator require a special pre-processing.
#' In practical terms, creating an expression of zero allows us to have better
#' control over the function call as suggested in the [rlang::exec()]
#' documentation. For instance, we can see how the function itself is being
#' called. Therefore, if an error occurs in one of the statistics, we will
#' have a direct traceback to the problematic call, as opposed to what happens
#'  directly using [rlang::exec()].
#'
#' @inheritParams decouple
#' @param fn Expression containing the name of the function to execute.
#' @param args Extra arguments to pass to the statistician under evaluation.
#'
#' @keywords internal
#' @noRd
.invoke_statistic <- function(fn,
                              args,
                              mat,
                              network,
                              .source,
                              .target,
                              mat_symbol,
                              network_symbol,
                              include_time,
                              minsize,
                              show_toy_call) {
    .toy_call <- expr(
        (!!fn)(
            mat = {{ mat_symbol }},
            network = {{ network_symbol }},
            .source = {{ .source }},
            .target = {{ .target }},
            !!!args)
    )

    if (show_toy_call) {
        utils::capture.output(rlang::qq_show(!!.toy_call)) %>%
            stringr::str_replace_all(pattern = "= \\^", "= ") %>%
            rlang::inform()
    }

    .call <- expr(
        (!!fn)(
            mat = mat,
            network = network,
            .source = {{ .source }},
            .target = {{ .target }},
            !!!args)
    )

    if (include_time) {
        .start_time <- Sys.time()
        eval(.call) %>%
            add_column(
                statistic_time = difftime(Sys.time(), .start_time),
                .after = "score"
            )
    } else {
        eval(.call)
    }
}

#' Convert object to symbol expression
#'
#' @param x An object or expression to convert to symbol
#'
#' @keywords internal
#' @noRd
.label_expr <- function(x) rlang::get_expr(enquo(x))
