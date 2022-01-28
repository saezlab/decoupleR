% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-dataset-converters.R
\name{extract_sets}
\alias{extract_sets}
\title{Extract sets}
\usage{
extract_sets(network)
}
\arguments{
\item{network}{Tibble or dataframe with edges and it's associated metadata.}
}
\description{
Extracts feature sets from a renamed network (see \link{rename_net}).
}
\examples{
inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
mat <- readRDS(file.path(inputs_dir, "input-expr_matrix.rds"))
network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
network <- rename_net(network, tf, target, mor, likelihood)
extract_sets(network)
}