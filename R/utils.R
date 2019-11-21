#' Tidy a matrix
#'
#' This utility function takes a matrix and converts it to a tidy format and
#' adds if available observations' meta data.
#'
#' @param mat A matrix with observations/features in rows and variables in
#' columns
#' @param feature Class name of observations/features, e.g.
#' transcription_factors
#' @param key Class name of variables, e.g. samples
#' @param value Class name of matrix values, e.g. activities
#' @param meta Data frame with meta data of the observations. To map the meta
#' data to the tidied table the observation/feature column name must be
#' identical.
#'
#' @return Tidy table.
#'
#' @export
tdy = function(mat, feature, key, value, meta = NULL) {
  mat %>%
    data.frame(check.names = F, stringsAsFactors = F) %>%
    rownames_to_column(feature) %>%
    as_tibble() %>%
    gather({{key}}, {{value}}, -{{feature}}) %>%
    left_join(meta, by=feature)
}

#' Untidy a tibble
#'
#' This utility function takes a tidy tibble and converts it to a matrix.
#'
#' @param tbl A tidy tibble
#' @param feature Class name of observations/features present in tidy tibble
#' @param key Class name of key present in tidy tibble
#' @param value Class name of values in tidy tibble
#'
#' @return Matrix with observation in rows and variables in columns.
#'
#' @export
untdy = function(tbl, feature, key, value) {
  tbl %>%
    select({{feature}}, {{key}}, {{value}}) %>%
    spread({{key}}, {{value}}) %>%
    data.frame(row.names = 1, check.names = F, stringsAsFactors = F)
}
