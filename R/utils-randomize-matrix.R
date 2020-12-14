#' Randomize matrix
#'
#' Utility function used in functions that require permutations of the
#' expression matrix
#'
#' @param mat Matrix to randomize.
#' @param randomize_type How to randomize.
#'
#' @return Randomized matrix
#' @export
#' @examples
#' mat <- matrix(seq_len(9), ncol = 3)
#' mat
#'
#' set.seed(42)
#' randomize_matrix(mat, randomize_type = "rows")
#'
#' set.seed(42)
#' randomize_matrix(mat, randomize_type = "cols_independently")
randomize_matrix <- function(
    mat,
    randomize_type = c("rows", "cols_independently")) {
    randomize_type <- match.arg(randomize_type)

    switch(randomize_type,
        rows = mat[sample(nrow(mat)), ],
        cols_independently = apply(mat, 2, sample)
    ) %>%
        `row.names<-`(rownames(mat))
}
