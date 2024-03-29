% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replaceAlignmentCharacter.R
\name{replaceAlignmentCharacter}
\alias{replaceAlignmentCharacter}
\title{replaceAlignmentCharacter}
\usage{
replaceAlignmentCharacter(
  alignment = NULL,
  char.find = "N",
  char.replace = "-"
)
}
\arguments{
\item{alignment}{alignment in DNAbin, DNAStringSet, list, and matrix formats}

\item{char.find}{the character to find}

\item{char.replace}{the character to replace "char.find" with}
}
\value{
returns DNAStringSet of alignment with replaced characters
}
\description{
Function for replacing a character with another in a multiple sequence alignment
}
\examples{

your.tree = ape::read.tree(file = "file-path-to-tree.tre")
astral.data = astralPlane(astral.tree = your.tree,
                          outgroups = c("species_one", "species_two"),
                          tip.length = 1)

}
