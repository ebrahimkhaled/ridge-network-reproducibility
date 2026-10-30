\name{mps_select}
\alias{mps_select}
\title{Multi-penalty permutation selection of edges}
\description{
Selects the pairs whose ridge partial correlation exceeds a permutation
threshold at some penalty in a grid. For each penalty the threshold is an
order statistic of the maximal absolute partial correlation over the
observed data and \code{B} data sets in which each column is permuted
independently, at level \code{alpha} divided by the number of penalties.
}
\usage{
mps_select(x, mult = c(0.3, 1, 3, 10), B = 99L, alpha = 0.1)
}
\arguments{
\item{x}{numeric data matrix with observations in rows.}
\item{mult}{penalty grid in units of \eqn{\mathrm{tr}(S) / p}.}
\item{B}{number of permutations.}
\item{alpha}{level in \eqn{(0, 1)}.}
}
\value{A list with the logical \code{adjacency} matrix, the penalties
\code{lambda} and the thresholds \code{threshold}.}
\details{If the rows of \code{x} are independent and identically distributed
and its columns are mutually independent, the probability that any edge is
selected is at most \code{alpha}, for every sample size and dimension. No
edge can be selected unless \code{B + 1 >= length(mult) / alpha}.}
\examples{
set.seed(1)
x <- matrix(rnorm(60 * 8), 60, 8)
x[, 2] <- x[, 1] + rnorm(60, sd = 0.5)
sel <- mps_select(x, B = 49)
which(sel$adjacency & upper.tri(sel$adjacency), arr.ind = TRUE)
}
