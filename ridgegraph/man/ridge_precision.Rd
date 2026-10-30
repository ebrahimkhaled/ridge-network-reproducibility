\name{ridge_precision}
\alias{ridge_precision}
\alias{ridge_pcor}
\title{Ridge precision matrix and ridge partial correlations}
\description{
Computes the ridge precision estimator \eqn{(S + \lambda I)^{-1}} from a data
matrix and the corresponding partial correlations
\eqn{-q_{ij} / (q_{ii} q_{jj})^{1/2}}.
}
\usage{
ridge_precision(x, lambda)
ridge_pcor(x, lambda)
}
\arguments{
\item{x}{numeric data matrix with observations in rows. Columns are centred
internally and \eqn{S = X^\top X / (n - 1)}.}
\item{lambda}{a single positive penalty.}
}
\value{A symmetric \eqn{p \times p} matrix. For \code{ridge_pcor} the
diagonal is set to one.}
\details{The partial correlations do not depend on the scale of the
estimator. For Gaussian data they estimate the partial correlations of
\eqn{(\Sigma + \kappa I)^{-1}}, where \eqn{\kappa} is the effective penalty
returned by \code{\link{effective_penalty}}.}
\examples{
set.seed(1)
x <- matrix(rnorm(40 * 10), 40, 10)
r <- ridge_pcor(x, lambda = 0.5)
range(r[upper.tri(r)])
}
