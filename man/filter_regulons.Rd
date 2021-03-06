% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/decoupleR-pre.R
\name{filter_regulons}
\alias{filter_regulons}
\title{Filter network by size of regulons}
\usage{
filter_regulons(network, .source, min_size = 1, max_size = Inf)
}
\arguments{
\item{network}{Tibble or dataframe with edges and it's associated metadata.}

\item{.source}{Column with source nodes.}

\item{min_size}{Minimum number of targets allowed per regulon.}

\item{max_size}{Maximum number of targets allowed per regulon.}
}
\value{
Filtered tibble.
}
\description{
Keep only sources which satisfied the condition \verb{min_size >= n <= max_size},
where \code{n} denotes the number of targets per source.
}
\examples{
inputs_dir <- system.file("testdata", "inputs", package = "decoupleR")
network <- readRDS(file.path(inputs_dir, "input-dorothea_genesets.rds"))
filter_regulons(network, .source = tf, min_size = 30, max_size = 50)
}
