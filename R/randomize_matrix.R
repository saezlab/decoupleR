#' Randomize matrix
#'
#' @param mat Matrix to randomize.
#' @param randomize_type How to randomize.
#'
#' @return Randomized matrix
#' @export
#'
randomize_matrix <- function(mat, randomize_type = c("rows", "cols_independently")) {
  randomize_type <- match.arg(randomize_type)

  switch(randomize_type,
    rows = mat[sample(nrow(mat)), ],
    cols_independently = apply(mat, 2, sample)
  ) %>%
    `row.names<-`(rownames(mat))
}
