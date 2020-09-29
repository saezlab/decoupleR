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
#' @keywords internal
tdy = function(mat, feature, key, value, meta = NULL) {
  mat %>%
    data.frame(check.names = FALSE, stringsAsFactors = FALSE) %>%
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
#' @keywords internal
untdy = function(tbl, feature, key, value) {
  tbl %>%
    select({{feature}}, {{key}}, {{value}}) %>%
    spread({{key}}, {{value}}) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
}




#' Helper function
#'
#' Helper function to convert between commonly used gene IDs
#' (before committing I should create a resource separate from annotation dbi)
#'
#' @param de_data Table with de results
#' @param input_type input gene IDs e.g. ENSEMBL
#' @param output_type output gene IDs e.g. Symbol
#' @export
#' @return Returns a df with gene names convert to symbols
convert_gene_type <- function(de_data, input_type, output_type) {
  geneIDs1 <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys=de_data$X,
                                    keytype = input_type,
                                    columns = c(output_type, input_type))

  geneIDs1 <- subset(geneIDs1, (!duplicated(geneIDs1[[output_type]])))
  geneIDs1 <- subset(geneIDs1, !is.na(geneIDs1[[output_type]]))
  data_me <- merge(de_data, geneIDs1, by.x = "X", by.y = input_type)
  data_me$X <- data_me[[output_type]]
  data_me <- within(data_me, rm(list=sub("[.]test","",output_type)))
  return(data_me)
}


