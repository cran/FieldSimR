% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field_trial_error.R
\name{field_trial_error}
\alias{field_trial_error}
\title{Simulate plot errors in plant breeding field trials}
\usage{
field_trial_error(
  ntraits = 1,
  nenvs = 1,
  nblocks = 2,
  block.dir = "col",
  ncols = 10,
  nrows = 20,
  varR = 1,
  ScorR = NULL,
  RcorR = NULL,
  EcorR = NULL,
  spatial.model = "Bivariate",
  complexity = NULL,
  plot.length = 8,
  plot.width = 2,
  col.cor = 0.5,
  row.cor = 0.7,
  prop.spatial = 0.5,
  ext.ord = "random",
  ext.dir = "row",
  prop.ext = 0,
  return.effects = FALSE
)
}
\arguments{
\item{ntraits}{Number of traits to be simulated.}

\item{nenvs}{Number of environments to be simulated.}

\item{nblocks}{A vector defining the number of blocks in each environment.
If only one value is specified, all environments will be assigned the same number.}

\item{block.dir}{A vector defining the block direction in each environment. Use
'col' for side-by-side (default), 'row' for above-and-below, or
NA if only one block is simulated. If only one value is specified, all environments
will be assigned the same block direction.}

\item{ncols}{A vector defining the number of columns in each environment. If only one
value is specified, all environments will be assigned the same number.}

\item{nrows}{A vector defining the number of rows in each environment. If only one
value is specified, all environments will be assigned the same number.}

\item{varR}{A vector of error variances for each environment-within-trait combination. If only
one value is specified, all combinations will be assigned the same error variance.}

\item{ScorR}{A matrix of spatial error correlations between traits. If not specified and
spatial trend is simulated, a diagonal matrix is constructed.}

\item{RcorR}{A matrix of random error correlations between traits. If not specified and
random error is simulated, a diagonal matrix is constructed.}

\item{EcorR}{A matrix of extraneous error correlations between traits. If not specified and
extraneous variation is simulated, a diagonal matrix is constructed. \cr
\strong{Note:} the same correlation between traits is used for the column and row errors.
Currently only implemented when \code{ext.ord = "random"}.}

\item{spatial.model}{A character string defining the model used to simulate spatial trend.
Use 'Bivariate' for bivariate interpolation (default) or 'AR1' for a separable first-order
autoregressive process. Bivariate interpolation is implemented with the \code{interp}
function of the R package \href{https://CRAN.R-project.org/package=interp}{`interp`}.}

\item{complexity}{A vector defining the complexity of the simulated spatial trend in each
environment when \code{spatial.model = "Bivariate"}. If only one value is specified,
all environments will be assigned the same complexity.
If not specified and \code{spatial.model = "Bivariate"}, the complexity is set to half the maximum
number of columns and rows in each environment.}

\item{plot.length}{A vector of plot lengths for each environment (column direction). If only one value is specified,
all environments will be assigned the same plot length. Only required when \code{spatial.model = "Bivariate"}.}

\item{plot.width}{A vector of plot widths for each environment (row direction). If only one value is specified,
all environments will be assigned the same plot width. Only required when \code{spatial.model = "Bivariate"}.}

\item{col.cor}{A vector of column autocorrelations for each environment.
If only one value is specified, all environments will be assigned the same column autocorrelation.
Only required when \code{spatial.model = "AR1"}.}

\item{row.cor}{A vector of row autocorrelations for each environment.
If only one value is specified, all environments will be assigned the same row
autocorrelation. Only required when \code{spatial.model = "AR1"}.}

\item{prop.spatial}{A vector defining the proportion of spatial trend for each environment-within-trait
combination. If only one value is specified, all combinations will be assigned the proportion.}

\item{ext.ord}{A character string defining the method used to simulate extraneous variation.
Use 'random' (default) for random variation between neighbouring columns and/or rows or
'zig-zag' for alternating positive and negative values.}

\item{ext.dir}{A vector defining the direction of extraneous variation for each environment.
Use 'row' (default) for row variation, 'col' for column variation, 'both' for variation in both directions,
or NA if no extraneous variation is simulated. When \code{ext.dir = "both"}, half the variance is assigned to the columns and
half is assigned to the rows. If only one value is specified, all environments will be
assigned the same direction.}

\item{prop.ext}{A vector defining the proportion of extraneous variation for each environment-within-trait
combination. If only one value is specified, all combinations will be assigned the same proportion.}

\item{return.effects}{When \code{TRUE} (default is \code{FALSE}), a list is returned with additional entries
containing the spatial, random, and extraneous error terms for each trait.}
}
\value{
A data frame with columns 'env', 'block', 'col', and 'row', followed by the
  simulated plot errors for each trait. When \code{return.effects = TRUE}, a list is returned with additional entries
  containing the spatial, random, and extraneous error terms for each trait.
}
\description{
Creates a data frame of simulated plot errors in multi-environment field trials for one or
more traits. The plot errors capture spatial trend, random error (noise), and extraneous variation.
Spatial trend is simulated using bivariate interpolation or a separable first-order
autoregressive (AR1) process. Random error is simulated using an independent process. Extraneous
variation is simulated using random or zig-zag ordering between neighbouring columns and/or rows.
The three error components are combined at a user-defined ratio. \cr
Correlated plot errors can be simulated between traits by setting different correlation structures
for each error component. A separable structure is assumed between traits and plots within environments,
but different error variances can be specified for each environment-within-trait combination.
}
\examples{
# Simulate plot errors for two traits in two environments using an AR1 model
# for spatial variation.

# Error variances for the four environment-within-trait combinations.
varR <- c(0.2, 0.4, 10, 15) # Trait 1 x 2 environments, Trait 2 x 2 environments

# Spatial error correlations between the two simulated traits.
ScorR <- matrix(c(
  1.0, 0.2,
  0.2, 1.0
), ncol = 2)

error_ls <- field_trial_error(
  ntraits = 2,
  nenvs = 2,
  nblocks = 2,
  block.dir = "row",
  ncols = 10,
  nrows = 20,
  varR = varR,
  ScorR = ScorR,
  spatial.model = "AR1",
  col.cor = 0.5,
  row.cor = 0.7,
  prop.spatial = 0.4,
  ext.ord = "zig-zag",
  ext.dir = "row",
  prop.ext = 0.2,
  return.effects = TRUE
)

}
