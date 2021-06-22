#' Tidy eval helpers
#'
#' @description
#'
#' * [rlang::sym()] creates a symbol from a string and
#'   [`syms()`][rlang::sym] creates a list of symbols from a
#'   character vector.
#'
#' * [`enquo()`][rlang::nse-defuse] and
#'   [`enquos()`][rlang::nse-defuse] delay the execution of one or
#'   several function arguments. `enquo()` returns a single quoted
#'   expression, which is like a blueprint for the delayed computation.
#'   `enquos()` returns a list of such quoted expressions.
#'
#' * [`expr()`][rlang::nse-defuse] quotes a new expression _locally_. It
#'   is mostly useful to build new expressions around arguments
#'   captured with [enquo()] or [enquos()]:
#'   `expr(mean(!!enquo(arg), na.rm = TRUE))`.
#'
#' * [rlang::as_name()] transforms a quoted variable name
#'   into a string. Supplying something else than a quoted variable
#'   name is an error.
#'
#'   That's unlike [rlang::as_label()] which also returns
#'   a single string but supports any kind of R object as input,
#'   including quoted function calls and vectors. Its purpose is to
#'   summarise that object into a single label. That label is often
#'   suitable as a default name.
#'
#'   If you don't know what a quoted expression contains (for instance
#'   expressions captured with `enquo()` could be a variable
#'   name, a call to a function, or an unquoted constant), then use
#'   `as_label()`. If you know you have quoted a simple variable
#'   name, or would like to enforce this, use `as_name()`.
#'
#' To learn more about tidy eval and how to use these tools, visit
#' <https://tidyeval.tidyverse.org> and the
#' [Metaprogramming section](https://adv-r.hadley.nz/metaprogramming.html) of
#' [Advanced R](https://adv-r.hadley.nz).
#'
#' @md
#' @name tidyeval
#' @keywords internal
#' @importFrom rlang expr enquo enquos sym syms .data := as_name as_label quo_is_null quo_is_missing abort exec
#' @aliases expr enquo enquos sym syms .data := as_name as_label quo_is_null quo_is_missing abort exec
#' @export expr enquo enquos sym syms .data := as_name as_label quo_is_null quo_is_missing abort exec
#' @examples
#' if (FALSE) {
#'     help("nse-defuse", package = "rlang")
#' }
NULL
