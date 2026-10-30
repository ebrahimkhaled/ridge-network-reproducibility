## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## -----------------------------------------------------------------------------
library(ridgegraph)
set.seed(1)
p <- 30; n <- 20
theta <- diag(p)
for (i in 1:(p - 1)) theta[i, i + 1] <- theta[i + 1, i] <- 0.4
x <- matrix(rnorm(n * p), n, p) %*% chol(solve(theta))
effective_penalty(x, c(0.05, 0.5, 5))
effective_penalty_floor(x)

## -----------------------------------------------------------------------------
sel <- mps_select(x, B = 99)
sum(sel$adjacency[upper.tri(sel$adjacency)])

