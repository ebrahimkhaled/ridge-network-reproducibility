###############################################################################
## Fixed-point effective penalty and Ky Fan floor for Table 1 of
##   "What does a ridge network estimate? Effective penalty and leakage in
##    ridge precision matrix estimation"
##
## Author : M. Arashi
## Purpose: adds two deterministic columns to the effective-penalty grid of
##          Table 1 (band model, p = 200):
##            kappa_bar  root of  lambda / k + n^{-1} sum_j s_j / (s_j + k) = 1,
##                       s_j the eigenvalues of Sigma (fixed-point lemma);
##            kyfan      lambda + n^{-1} times the sum of the (p - n + 1)_+
##                       smallest eigenvalues of Sigma (floor of the effective
##                       penalty).
##          Nothing here is random; the Monte Carlo columns are read from
##          sim_kappa_grid.csv, written by sim_ridge_support.R.
## Requirements: R >= 4.1; package huge. An extra library folder can be given
##               in the environment variable RIDGE_RLIB.
## Usage  : Rscript table1_fixed_point.R   (run from code/)
## Output : ../results/sim_kappa_fixedpoint.csv
###############################################################################

lib <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(huge))
OUT <- ".."

## band model of Table 1 (the huge band generator uses no random numbers)
p <- 200
Sigma <- as.matrix(huge.generator(d = p, graph = "band", g = 1, v = 0.4, u = 0.1,
                                  verbose = FALSE)$sigma)
ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values

kappa_bar <- function(n, lambda) {
  psi <- function(k) lambda / k + sum(ev / (ev + k)) / n - 1
  uniroot(psi, c(lambda, lambda + sum(ev) / n + 1e-12), tol = 1e-14)$root
}
kyfan_floor <- function(n, lambda) {
  k <- max(0, p - n + 1)
  lambda + if (k > 0) sum(sort(ev)[1:k]) / n else 0
}

tab <- read.csv(file.path(OUT, "results", "sim_kappa_grid.csv"))
tab$kappa_bar <- mapply(kappa_bar, tab$n, tab$lambda)
tab$kyfan <- mapply(kyfan_floor, tab$n, tab$lambda)
tab$relerr_bar <- (tab$kappa_bar - tab$kappa_mc) / tab$kappa_mc
write.csv(tab, file.path(OUT, "results", "sim_kappa_fixedpoint.csv"), row.names = FALSE)

cat(sprintf("largest |kappa_bar - kappa_mc| / kappa_mc: %.4f\n", max(abs(tab$relerr_bar))))
cat(sprintf("Ky Fan floor below the Monte Carlo value in every cell: %s\n",
            all(tab$kyfan <= tab$kappa_mc)))
