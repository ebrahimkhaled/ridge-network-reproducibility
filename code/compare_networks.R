## Comparing ridge networks from samples of different size: at a common
## nominal penalty or at a common estimated effective penalty?
##
## Section 6.2 of the article: quarterly against yearly networks of the stock
## data, and the semi-synthetic designs with a known ordering.
##
## The semi-synthetic designs weaken the correlations BETWEEN sectors of one
## population.  We do not use Sigma_B = (1-w) Sigma_A + w I: that matrix has
## exactly the partial correlations of A at a larger effective penalty, so a
## common effective penalty would be right by construction.
##
## Decision rule, fixed before the run: the claim "a common nominal penalty
## misorders networks from samples of different size, and a common effective
## penalty removes most of that error" stands only if, in every design, the
## nominal comparison gives the population ordering in at most 50% of
## replicates and the matched comparison in at least 80% at the largest
## difference.
##
## Strength of a network = mean signed partial correlation over all pairs.
##
## Output : ../results/compare_*.csv, sessionInfo_compare.txt

lib <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib)) .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(huge))
out <- file.path("..", "results"); dir.create(out, showWarnings = FALSE)

spec_S <- function(S) { e <- eigen(S, symmetric = TRUE); list(vec = e$vectors, val = pmax(e$values, 0)) }
ridge_Q <- function(sp, lambda) sp$vec %*% (t(sp$vec) / (sp$val + lambda))
pcor_of <- function(A) { d <- 1 / sqrt(diag(A)); R <- -A * outer(d, d); diag(R) <- 1; R }
kappa_hat <- function(sp, lambda, n) lambda / (1 - sum(sp$val / (sp$val + lambda)) / n)
lambda_for_kappa <- function(sp, n, k) uniroot(function(l) kappa_hat(sp, l, n) - k, c(1e-10, 1e6), tol = 1e-12)$root
msig <- function(R) mean(R[upper.tri(R)])
mabs <- function(R) mean(abs(R[upper.tri(R)]))

data(stockdata, package = "huge")
R0 <- diff(log(stockdata$data)); R0[abs(R0) > 0.25] <- 0
sector <- stockdata$info[, 2]
normal_scores <- function(Y) {
  Z <- apply(Y, 2, function(x) qnorm(rank(x, ties.method = "average") / (length(x) + 1)))
  scale(Z, scale = FALSE)
}

## ---------------------------------------------------------------- real data
fitw <- function(idx) {
  X <- normal_scores(R0[idx, ]); n <- nrow(X) - 1; S <- crossprod(X) / n
  list(sp = spec_S(S), n = n, corr = msig(cov2cor(S)), idx = idx)
}
years <- lapply(1:5, function(w) fitw(((w - 1) * 251 + 1):(w * 251)))
quarters <- lapply(1:19, function(q) fitw(((q - 1) * 63 + 1):(q * 63)))
real <- list()
for (LAM in c(0.4, 1)) for (q in seq_along(quarters)) for (y in seq_along(years)) {
  Q <- quarters[[q]]; Y <- years[[y]]
  kq <- kappa_hat(Q$sp, LAM, Q$n); ly <- lambda_for_kappa(Y$sp, Y$n, kq)
  Rq <- pcor_of(ridge_Q(Q$sp, LAM)); Ryn <- pcor_of(ridge_Q(Y$sp, LAM)); Rym <- pcor_of(ridge_Q(Y$sp, ly))
  real[[length(real) + 1]] <- data.frame(lambda = LAM, quarter = q, year = y,
    overlap = length(intersect(Q$idx, Y$idx)) > 0, kappa_q = kq, kappa_y = kappa_hat(Y$sp, LAM, Y$n),
    nominal = msig(Rq) > msig(Ryn), matched = msig(Rq) > msig(Rym), corr = Q$corr > Y$corr,
    abs_nominal = mabs(Rq) > mabs(Ryn), abs_matched = mabs(Rq) > mabs(Rym))
}
real <- do.call(rbind, real)
write.csv(real, file.path(out, "compare_real.csv"), row.names = FALSE)
for (LAM in c(0.4, 1)) for (ov in c(FALSE, TRUE)) {
  d <- real[real$lambda == LAM & (ov | !real$overlap), ]
  cat(sprintf("REAL lambda=%.1f %s pairs=%d: quarter stronger nominal %d matched %d corr %d | reversals %d | agree corr nominal %.2f matched %.2f | mean|rho| quarter stronger nominal %d matched %d\n",
      LAM, if (ov) "all" else "no-overlap", nrow(d), sum(d$nominal), sum(d$matched), sum(d$corr),
      sum(d$nominal != d$matched), mean(d$nominal == d$corr), mean(d$matched == d$corr),
      sum(d$abs_nominal), sum(d$abs_matched)))
}
## reversals per year (clustering), lambda 0.4, all pairs
d <- real[real$lambda == 0.4, ]; print(tapply(d$nominal != d$matched, d$year, sum))

## ----------------------------------------------------------- semi-synthetic
## The sample correlation matrices of 628 rows have tiny eigenvalues, and
## weakening blocks of them breaks positive definiteness.  Both halves are
## therefore lightly regularized in the same way, (R + 0.1 I) / 1.1, before
## use; this only defines the populations.
reg <- function(R) (R + 0.1 * diag(ncol(R))) / 1.1
X1 <- normal_scores(R0[1:628, ]); Sig1 <- reg(cov2cor(crossprod(X1) / (nrow(X1) - 1)))
X2 <- normal_scores(R0[629:1256, ]); Sig2 <- reg(cov2cor(crossprod(X2) / (nrow(X2) - 1)))
p <- ncol(Sig1)
pop <- function(Sig, k) msig(pcor_of(solve(Sig + k * diag(p))))
kgrid <- c(0, 0.5, 1, 2, 4, 8)
same <- outer(sector, sector, "==")
## B = (1 - w) A + w blockdiag(A): between-sector correlations times (1 - w),
## within-sector ones unchanged; a convex combination of positive definite
## matrices, so positive definite, and not a shift of the effective penalty
weaken <- function(Sig, w) { B <- Sig; B[!same] <- (1 - w) * B[!same]; B }

designs <- list(
  list(name = "between-sector correlations x0.9", A = Sig1, B = weaken(Sig1, 0.1)),
  list(name = "between-sector correlations x0.8", A = Sig1, B = weaken(Sig1, 0.2)),
  list(name = "between-sector correlations x0.7", A = Sig1, B = weaken(Sig1, 0.3)),
  list(name = "two halves of the data", A = Sig1, B = Sig2))
LAM <- 0.4; RREP <- 200; set.seed(20261125)
semi <- list(); reps <- list()
for (D in designs) {
  stopifnot(min(eigen(D$B, symmetric = TRUE, only.values = TRUE)$values) > 0)
  pd <- sapply(kgrid, function(k) pop(D$A, k) - pop(D$B, k))
  cat(sprintf("%s: population strength A - B at kappa %s = %s\n", D$name,
      paste(kgrid, collapse = "/"), paste(signif(pd, 3), collapse = " ")))
  ## the stronger population is observed with the SHORT sample, so the bias of
  ## a common nominal penalty works against the truth
  if (all(pd < 0)) { tmp <- D$A; D$A <- D$B; D$B <- tmp; pd <- -pd }
  if (!(all(pd > 0))) { cat("  ordering depends on kappa: design skipped\n"); next }
  LA <- chol(D$A); LB <- chol(D$B)
  for (r in 1:RREP) {
    xa <- scale(matrix(rnorm(63 * p), 63) %*% LA, scale = FALSE)
    xb <- scale(matrix(rnorm(251 * p), 251) %*% LB, scale = FALSE)
    spa <- spec_S(crossprod(xa) / 62); spb <- spec_S(crossprod(xb) / 250)
    ka <- kappa_hat(spa, LAM, 62)
    sa <- msig(pcor_of(ridge_Q(spa, LAM)))
    sbn <- msig(pcor_of(ridge_Q(spb, LAM)))
    sbm <- msig(pcor_of(ridge_Q(spb, lambda_for_kappa(spb, 250, ka))))
    reps[[length(reps) + 1]] <- data.frame(design = D$name, rep = r, kappa_A = ka,
      kappa_B = kappa_hat(spb, LAM, 250), diff_nominal = sa - sbn, diff_matched = sa - sbm)
  }
  rr <- do.call(rbind, reps[sapply(reps, function(z) z$design[1] == D$name)])
  popk <- pop(D$A, mean(rr$kappa_A)) - pop(D$B, mean(rr$kappa_A))   # population difference at the matched level
  semi[[length(semi) + 1]] <- data.frame(design = D$name, kappa_A = mean(rr$kappa_A), kappa_B = mean(rr$kappa_B),
    pop_diff_at_kappaA = popk, mean_diff_nominal = mean(rr$diff_nominal), mean_diff_matched = mean(rr$diff_matched),
    sd_diff_matched = sd(rr$diff_matched),
    correct_nominal = mean(rr$diff_nominal > 0), correct_matched = mean(rr$diff_matched > 0),
    mcse_matched = sqrt(mean(rr$diff_matched > 0) * (1 - mean(rr$diff_matched > 0)) / RREP))
  print(tail(semi, 1)[[1]], digits = 3)
  write.csv(do.call(rbind, reps), file.path(out, "compare_semisynthetic_replicates.csv"), row.names = FALSE)
}
write.csv(do.call(rbind, semi), file.path(out, "compare_semisynthetic.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_compare.txt"))
