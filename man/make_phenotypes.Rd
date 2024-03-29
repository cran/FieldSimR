% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_phenotypes.R
\name{make_phenotypes}
\alias{make_phenotypes}
\title{Simulate phenotypes - Combine simulated genetic values and plot errors}
\usage{
make_phenotypes(gv_df, error_df, randomise = FALSE)
}
\arguments{
\item{gv_df}{A data frame of genetic values. Must contain the columns 'env', 'rep', and 'id',
followed by the genetic values for each trait.}

\item{error_df}{A data frame of plot errors. Must contain the columns 'env', 'block',
'col', and 'row', followed by the plot errors for each trait.}

\item{randomise}{When TRUE, genotypes are randomly allocated to plots within blocks to generate
a randomized complete block design (RCBD).\cr
\strong{Note:} Other experimental designs must be generated externally.}
}
\value{
A data frame with columns 'env', 'block', 'column', 'row' and 'genotype', followed by
  the phenotypes for each trait.
}
\description{
Creates a data frame of simulated phenotypes for one or more traits by
combining simulated plot errors with genetic values (e.g. true, simulated or predicted).
The genetic values can be generated externally, but note that they must be stored in a
data frame as described below.
}
\examples{
# Simulate phenotypes by combining the genetic values and plot errors provided in
# the two example data frames 'df_gv_unstr' and 'df_error_bivar'.

gv_df <- df_gv_unstr
error_df <- df_error_bivar

pheno_df <- make_phenotypes(
  gv_df,
  error_df,
  randomise = TRUE
)
}
