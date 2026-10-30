\name{effective_penalty}
\alias{effective_penalty}
\alias{effective_penalty_floor}
\alias{penalty_for_kappa}
\title{Effective penalty of the ridge precision estimator}
\description{
Estimates the effective penalty \eqn{\kappa = \lambda (1 + \delta)} of the
ridge precision estimator through the degrees-of-freedom ratio
\eqn{\hat r = \mathrm{tr}\{S (S + \lambda I)^{-1}\} / (n - 1)}, as
\eqn{\hat\kappa = \lambda / (1 - \hat r)}; returns the limit of
\eqn{\hat\kappa} as \eqn{\lambda \to 0}; and finds the penalty that attains a
prescribed effective penalty. The columns of \code{x} are centred and
\eqn{S = X^\top X / (n - 1)}, so the effective sample size in all formulas
is \eqn{n - 1}, where \eqn{n} is the number of rows.
}
\usage{
effective_penalty(x, lambda)
effective_penalty_floor(x, tol = 1e-10)
penalty_for_kappa(x, kappa)
}
\arguments{
\item{x}{numeric data matrix with observations in rows.}
\item{lambda}{vector of positive penalties.}
\item{tol}{relative tolerance below which an eigenvalue of \eqn{S} is
treated as zero.}
\item{kappa}{target effective penalty, larger than
\code{effective_penalty_floor(x)}.}
}
\value{
\code{effective_penalty} returns a data frame with columns \code{lambda},
\code{df_ratio} (\eqn{\hat r}), \code{scale} (\eqn{1 - \hat r}, an estimate
of \eqn{1 / (1 + \delta)}) and \code{kappa}. \code{effective_penalty_floor}
returns the harmonic mean of the positive eigenvalues of \eqn{S} when their
number equals the effective sample size, and zero otherwise.
\code{penalty_for_kappa} returns a penalty.
}
\details{The estimated effective penalty is continuous and strictly
increasing in \eqn{\lambda}, so \code{penalty_for_kappa} has a unique
solution.}
\examples{
set.seed(1)
x <- matrix(rnorm(30 * 60), 30, 60)
effective_penalty(x, c(0.01, 0.1, 1))
effective_penalty_floor(x)
}
