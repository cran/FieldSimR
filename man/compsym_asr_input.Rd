% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compound_symmetry_gxe.R
\name{compsym_asr_input}
\alias{compsym_asr_input}
\title{Simulate genetic values based on a compound symmetry model for GxE interaction - `AlphaSimR`
input parameters}
\usage{
compsym_asr_input(
  ntraits = 1,
  nenvs = 2,
  mean = 0,
  var = 1,
  prop.main = 0.5,
  corA = NULL,
  meanDD = NULL,
  varDD = NULL,
  prop.mainDD = NULL,
  corDD = NULL,
  relAA = NULL,
  prop.mainAA = NULL,
  corAA = NULL
)
}
\arguments{
\item{ntraits}{Number of traits to be simulated.}

\item{nenvs}{Number of environments to be simulated (minimum of two).}

\item{mean}{A vector of mean genetic values for each environment-within-trait combination.
If only one value is specified, all combinations will be assigned the same mean.}

\item{var}{A vector of genetic variances for each trait. \cr
\strong{Note:} When \code{useVarA = TRUE} is specified in `AlphaSimR` (default), the values in
\code{var} represent the additive genetic variances, otherwise they represent the
total (additive + non-additive) genetic variances.}

\item{prop.main}{A vector defining the proportion of main effect variance for each trait.
If only one value is specified, all traits will be assigned the same proportion. \cr
\strong{Note:} \code{0 < prop.main < 1}.}

\item{corA}{A matrix of additive genetic correlations between traits. By default, a diagonal
matrix is constructed.}

\item{meanDD}{A vector of mean dominance degrees for each environment-within-trait combination
(similar to \code{mean}). If only one value is specified, all combinations will be assigned
the same mean. By default, \code{meanDD = NULL} and dominance is not simulated.}

\item{varDD}{A vector of dominance degree variances for each trait.}

\item{prop.mainDD}{A vector defining the proportion of dominance degree main effect
variance for each trait (similar to \code{prop.main}).
If only one value is specified, all traits will be assigned the same proportion. \cr
\strong{Note:} \code{0 < prop.mainDD < 1}.}

\item{corDD}{A matrix of dominance degree correlations between traits (similar
to \code{corA}). If not specified and dominance is simulated, a diagonal matrix is constructed.}

\item{relAA}{A vector defining the relative magnitude of additive-by-additive (epistatic) variance
to additive genetic variance for each trait, that is in a diploid organism with
allele frequency of 0.5.
If only one value is specified, all traits will be assigned the same relative magnitude.}

\item{prop.mainAA}{A vector defining the proportion of epistatic main effect variance for each
trait (similar to \code{prop.main}). If only one value is specified, all traits will be assigned the
same proportion. \cr
\strong{Note:} \code{0 < prop.mainAA < 1}.}

\item{corAA}{A matrix of epistatic correlations between traits (similar to
\code{corA}). If not specified and epistasis is simulated, a diagonal matrix is constructed.}
}
\value{
A list with input parameters for `AlphaSimR`, which are used to simulate
  correlated genetic values based on a compound symmetry model for GxE interaction.
}
\description{
Creates a list of input parameters for
\href{https://CRAN.R-project.org/package=AlphaSimR}{`AlphaSimR`} to simulate
genetic values in multiple environments for one or more traits based on a compound symmetry
model for genotype-by-environment (GxE) interaction. \cr
This function utilises the ability of `AlphaSimR` to simulate correlated traits.
The wrapper function \code{compsym_asr_input()} is used to specify the input parameters required in `AlphaSimR`.
After simulating the genetic values, the wrapper function \link[FieldSimR]{compsym_asr_output} can be used to
generate a data frame with output values.
}
\details{
The compound symmetry model assumes the same genetic variance for each environment
and the same genetic covariance between each pair of environments. New functionality
is being implemented which relaxes the former assumption
(also see \link[FieldSimR]{unstr_asr_output}).

\strong{Note:} `AlphaSimR` can simulate different biological effects (see:
\href{https://gaynorr.github.io/AlphaSimR/reference/SimParam.html}{SimParam}).
\itemize{
  \item For additive traits use \code{addTraitA()}.
  \item For additive + dominance traits use \code{addTraitAD()}.
  \item For additive + epistatic traits use \code{addTraitAE()}.
  \item For additive + dominance + epistatic traits use \code{addTraitADE()}.
  }
Check the \code{useVarA} argument of these functions when simulating non-additive traits.
}
\examples{
# Simulate genetic values with 'AlphaSimR' for two additive + dominance traits
# in two environments based on a compound symmetry model.

# 1. Define the genetic architecture of the simulated traits.
# Mean genetic values and mean dominance degrees.
mean <- c(4.9, 5.4, 235.2, 228.5) # Trait 1 x 2 environments, Trait 2 x 2 environments
meanDD <- c(0.4, 0.4, 0.1, 0.1) # Trait 1 and 2, same value for both environments

# Additive genetic variances and dominance degree variances.
var <- c(0.08, 13) # Different values for Traits 1 and 2
varDD <- 0.2 # Same value for Traits 1 and 2

# Proportion of additive and dominance degree main effect variances.
prop.main <- c(0.4, 0.6) # Different values for Traits 1 and 2
prop.mainDD <- 0.4 # Same value for Traits 1 and 2

# Additive and dominance degree correlations between the two simulated traits.
corA <- matrix(c(
  1.0, 0.5,
  0.5, 1.0
), ncol = 2)
corDD <- diag(2) # Assuming independence

input_asr <- compsym_asr_input(
  ntraits = 2,
  nenvs = 2,
  mean = mean,
  var = var,
  prop.main = prop.main,
  corA = corA,
  meanDD = meanDD,
  varDD = varDD,
  prop.mainDD = prop.mainDD,
  corDD = corDD
)

}
