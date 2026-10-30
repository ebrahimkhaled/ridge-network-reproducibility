# Core functions of the ridgegraph package.
#
# The sample covariance is taken as S = X'X / n_eff, where X is the data
# matrix after column centring and n_eff = nrow(X) - 1. After centring, S is
# distributed as the sample covariance of n_eff independent Gaussian vectors,
# so n_eff is the sample size that enters the effective-penalty formulas.

.check_x <- function(x) {
  if (!is.matrix(x) || !is.numeric(x))
    stop("'x' must be a numeric matrix")
  if (anyNA(x)) stop("'x' must not contain missing values")
  if (nrow(x) < 3L || ncol(x) < 2L)
    stop("'x' must have at least 3 rows and 2 columns")
  invisible(TRUE)
}

.spec <- function(x) {
  xc <- scale(x, center = TRUE, scale = FALSE)
  neff <- nrow(x) - 1L
  e <- eigen(crossprod(xc) / neff, symmetric = TRUE)
  list(vec = e$vectors, val = pmax(e$values, 0), neff = neff, xc = xc)
}

.q <- function(sp, lambda) sp$vec %*% (t(sp$vec) / (sp$val + lambda))

.pcor <- function(a) {
  d <- 1 / sqrt(diag(a))
  r <- -a * outer(d, d)
  diag(r) <- 1
  r
}

# 1 - rhat written as (n - r) / n + (lambda / n) sum_{s_j > 0} 1 / (s_j + lambda),
# with r the number of positive eigenvalues, to avoid cancellation for small
# lambda
.shat <- function(sp, lambda, tol = 1e-10) {
  pos <- sp$val[sp$val > tol * max(sp$val)]
  (sp$neff - length(pos)) / sp$neff + lambda * sum(1 / (pos + lambda)) / sp$neff
}

.kappa <- function(sp, lambda) lambda / .shat(sp, lambda)

ridge_precision <- function(x, lambda) {
  .check_x(x)
  if (!is.numeric(lambda) || length(lambda) != 1L || lambda <= 0)
    stop("'lambda' must be a single positive number")
  .q(.spec(x), lambda)
}

ridge_pcor <- function(x, lambda) {
  .pcor(ridge_precision(x, lambda))
}

effective_penalty <- function(x, lambda) {
  .check_x(x)
  if (!is.numeric(lambda) || any(lambda <= 0))
    stop("'lambda' must contain positive numbers")
  sp <- .spec(x)
  sc <- vapply(lambda, function(l) .shat(sp, l), numeric(1))
  data.frame(lambda = lambda, df_ratio = 1 - sc, scale = sc,
             kappa = lambda / sc)
}

effective_penalty_floor <- function(x, tol = 1e-10) {
  .check_x(x)
  sp <- .spec(x)
  pos <- sp$val[sp$val > tol * max(sp$val)]
  if (length(pos) < sp$neff) 0 else sp$neff / sum(1 / pos)
}

penalty_for_kappa <- function(x, kappa) {
  .check_x(x)
  sp <- .spec(x)
  fl <- effective_penalty_floor(x)
  if (kappa <= fl)
    stop("'kappa' must exceed the limit of the effective penalty, ",
         format(fl, digits = 4))
  uniroot(function(l) .kappa(sp, l) - kappa, c(1e-12, 1e8),
          tol = 1e-12)$root
}

extrapolated_ridge <- function(x, lambda, ratio = 2) {
  .check_x(x)
  if (ratio <= 1) stop("'ratio' must exceed 1")
  sp <- .spec(x)
  k1 <- .kappa(sp, lambda)
  l2 <- uniroot(function(l) .kappa(sp, l) - ratio * k1,
                c(lambda, 1e8), tol = 1e-12)$root
  k2 <- .kappa(sp, l2)
  s1 <- .shat(sp, lambda)
  s2 <- .shat(sp, l2)
  w1 <- k2 / (k2 - k1)
  w2 <- -k1 / (k2 - k1)
  est <- w1 * s1 * .q(sp, lambda) + w2 * s2 * .q(sp, l2)
  structure(list(precision = est, pcor = .pcor(est),
                 lambda = c(lambda, l2), kappa = c(k1, k2),
                 weights = c(w1, w2)), class = "extrapolated_ridge")
}

mps_select <- function(x, mult = c(0.3, 1, 3, 10), B = 99L, alpha = 0.1) {
  .check_x(x)
  if (alpha <= 0 || alpha >= 1) stop("'alpha' must lie in (0, 1)")
  p <- ncol(x)
  ut <- upper.tri(diag(p))
  sp <- .spec(x)
  lams <- mult * sum(sp$val) / p
  robs <- lapply(lams, function(l) .pcor(.q(sp, l)))
  tobs <- vapply(robs, function(r) max(abs(r[ut])), numeric(1))
  tb <- matrix(0, B, length(lams))
  for (b in seq_len(B)) {
    spb <- .spec(apply(x, 2L, sample))
    tb[b, ] <- vapply(lams, function(l) max(abs(.pcor(.q(spb, l))[ut])),
                      numeric(1))
  }
  if (B + 1 < length(lams) / alpha)
    warning("no edge can be selected: B + 1 must be at least ",
            "length(mult) / alpha = ", ceiling(length(lams) / alpha))
  k <- ceiling((1 - alpha / length(lams)) * (B + 1))
  tau <- vapply(seq_along(lams), function(g) sort(c(tobs[g], tb[, g]))[k],
                numeric(1))
  adj <- matrix(FALSE, p, p)
  for (g in seq_along(lams)) adj <- adj | (abs(robs[[g]]) > tau[g])
  diag(adj) <- FALSE
  list(adjacency = adj, lambda = lams, threshold = tau)
}
