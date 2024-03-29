% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trimAlignmentColumns.R
\name{trimAlignmentColumns}
\alias{trimAlignmentColumns}
\title{trimAlignmentColumns}
\usage{
trimAlignmentColumns(alignment = NULL, min.gap.percent = 100)
}
\arguments{
\item{alignment}{alignment in DNAbin, DNAStringSet, list, and matrix formats}

\item{min.gap.percent}{minimum threshold gap percentage allowed to trim column}
}
\value{
returns DNAStringSet of column trimmed alignment
}
\description{
Function for trimming out alignment columns with too many gaps
}
\examples{

your.tree = ape::read.tree(file = "file-path-to-tree.tre")
astral.data = astralPlane(astral.tree = your.tree,
                          outgroups = c("species_one", "species_two"),
                          tip.length = 1)

}
