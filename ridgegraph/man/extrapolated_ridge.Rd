\name{extrapolated_ridge}
\alias{extrapolated_ridge}
\title{Extrapolation of ridge precision estimates across penalties}
\description{
Combines ridge estimates at two effective penalties \eqn{\kappa_1} and
\eqn{\kappa_2 = r \kappa_1} with weights \eqn{\kappa_2 / (\kappa_2 -
\kappa_1)} and \eqn{-\kappa_1 / (\kappa_2 - \kappa_1)}, after rescaling each
estimate by its estimated scale. The population version of this combination
differs from the precision matrix by a term of second order in the
penalties.
}
\usage{
extrapolated_ridge(x, lambda, ratio = 2)
}
\arguments{
\item{x}{numeric data matrix with observations in rows.}
\item{lambda}{first (smaller) penalty.}
\item{ratio}{ratio \eqn{\kappa_2 / \kappa_1 > 1} of the effective
penalties.}
}
\value{A list with the combined matrix \code{precision}, its partial
correlations \code{pcor}, the two penalties \code{lambda}, the two estimated
effective penalties \code{kappa} and the \code{weights}.}
\details{The sum of the absolute weights is \eqn{(r + 1) / (r - 1)}, which
multiplies the sampling variability.}
\examples{
set.seed(1)
x <- matrix(rnorm(50 * 20), 50, 20)
fit <- extrapolated_ridge(x, lambda = 0.5)
fit$kappa
}
