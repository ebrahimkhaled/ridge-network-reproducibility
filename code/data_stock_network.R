###############################################################################
## Real-data analysis for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Data   : `stockdata` from the R package huge (Zhao et al., 2012): daily
##          closing prices of 452 S&P 500 constituents over 1258 trading days
##          (2003-2008) and their GICS sector labels (10 sectors).
##          Source: huge::stockdata, available with the CRAN package huge.
## Purpose: (i) profile and clean the returns; (ii) estimate the effective
##          penalty of the ridge precision estimator in yearly windows
##          (p > n) and in the full sample (p < n); (iii) compare edge
##          selection layers on the same ridge fit and with the graphical
##          lasso and neighbourhood selection, through sector concordance and
##          window-to-window reproducibility; (iv) check the prediction that
##          ridge partial correlations approach marginal correlations as the
##          effective penalty grows.
## Requirements: R >= 4.1; packages huge, glmnet, fdrtool, igraph.
## Usage  : Rscript data_stock_network.R
## Output : ../results/data_*.csv and ../figs/fig_data_*.pdf
###############################################################################

suppressPackageStartupMessages({
  library(huge)
  library(glmnet)
  library(fdrtool)
  library(igraph)
})
OUT <- ".."
dir.create(file.path(OUT, "results"), showWarnings = FALSE)
dir.create(file.path(OUT, "figs"), showWarnings = FALSE)
set.seed(20260930)

NPERM <- 999                 # permutations for the calibrated threshold
ALPHA <- 0.1                 # family-wise level under complete independence
MULT  <- c(0.3, 1, 3, 10)    # penalty grid in units of tr(S) / p

## ---------------------------------------------------------------------------
## 1. Ridge building blocks (identical to the simulation script)
## ---------------------------------------------------------------------------
spec_S <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  list(vec = e$vectors, val = pmax(e$values, 0))
}
ridge_Q <- function(sp, lambda) sp$vec %*% (t(sp$vec) / (sp$val + lambda))
pcor_of <- function(A) {
  d <- 1 / sqrt(diag(A)); R <- -A * outer(d, d); diag(R) <- 1; R
}
## degrees-of-freedom ratio, scale and effective penalty; n is the effective
## sample size (number of rows minus one after centring)
r_hat <- function(sp, lambda, n) sum(sp$val / (sp$val + lambda)) / n
kappa_hat <- function(sp, lambda, n) lambda / (1 - r_hat(sp, lambda, n))
## limit of kappa_hat as lambda -> 0: the harmonic mean of the positive
## eigenvalues when rank(S) equals n, and 0 otherwise
kappa_floor_hat <- function(sp, n, tol = 1e-10) {
  pos <- sp$val[sp$val > tol * max(sp$val)]
  if (length(pos) < n) 0 else n / sum(1 / pos)
}
lambda_for_kappa <- function(sp, n, kappa_target) {
  uniroot(function(l) kappa_hat(sp, l, n) - kappa_target,
          c(1e-10, 1e6), tol = 1e-12)$root
}
richardson_Q <- function(sp, n, lambda1, cfac = 2) {
  k1 <- kappa_hat(sp, lambda1, n)
  lambda2 <- lambda_for_kappa(sp, n, cfac * k1)
  k2 <- kappa_hat(sp, lambda2, n)
  s1 <- 1 - r_hat(sp, lambda1, n); s2 <- 1 - r_hat(sp, lambda2, n)
  (k2 / (k2 - k1)) * s1 * ridge_Q(sp, lambda1) -
    (k1 / (k2 - k1)) * s2 * ridge_Q(sp, lambda2)
}
altridge_Q <- function(sp, lambda) {
  d <- sqrt(lambda + sp$val^2 / 4) + sp$val / 2
  sp$vec %*% (t(sp$vec) / d)
}
cv_lambda <- function(X, grid, K = 5, type = c("ridge", "alt")) {
  type <- match.arg(type)
  n <- nrow(X); fold <- sample(rep_len(1:K, n))
  loss <- matrix(0, K, length(grid))
  for (k in 1:K) {
    Xtr <- X[fold != k, , drop = FALSE]; Xte <- X[fold == k, , drop = FALSE]
    sp <- spec_S(crossprod(Xtr) / nrow(Xtr)); Ste <- crossprod(Xte) / nrow(Xte)
    sd_ <- colSums(sp$vec * (Ste %*% sp$vec))
    for (g in seq_along(grid)) {
      dv <- if (type == "ridge") 1 / (sp$val + grid[g])
            else 1 / (sqrt(grid[g] + sp$val^2 / 4) + sp$val / 2)
      loss[k, g] <- sum(sd_ * dv) - sum(log(dv))
    }
  }
  grid[which.min(colMeans(loss))]
}
lfdr_select <- function(R, cut = 0.9) {
  ut <- upper.tri(R)
  fit <- fdrtool(R[ut], statistic = "correlation", plot = FALSE,
                 verbose = FALSE)
  A <- matrix(FALSE, nrow(R), ncol(R)); A[ut] <- (1 - fit$lfdr) >= cut
  A | t(A)
}
mps_select <- function(X, mult = MULT, B = NPERM, alpha = ALPHA) {
  n <- nrow(X) - 1; p <- ncol(X); ut <- upper.tri(diag(p))
  S <- crossprod(X) / n
  lams <- mult * sum(diag(S)) / p
  sp <- spec_S(S)
  Robs <- lapply(lams, function(l) pcor_of(ridge_Q(sp, l)))
  Tobs <- sapply(Robs, function(R) max(abs(R[ut])))
  Tb <- matrix(0, B, length(lams))
  for (b in 1:B) {
    Xb <- scale(apply(X, 2, sample), center = TRUE, scale = FALSE)
    spb <- spec_S(crossprod(Xb) / n)
    Tb[b, ] <- sapply(lams, function(l) max(abs(pcor_of(ridge_Q(spb, l))[ut])))
  }
  k <- ceiling((1 - alpha / length(lams)) * (B + 1))
  A <- matrix(FALSE, p, p)
  tau <- sapply(seq_along(lams), function(g) sort(c(Tobs[g], Tb[, g]))[k])
  for (g in seq_along(lams)) A <- A | (abs(Robs[[g]]) > tau[g])
  diag(A) <- FALSE
  list(A = A, lambda = lams, tau = tau)
}
screen_refit <- function(X, R, dn) {
  p <- ncol(X); n <- nrow(X); B <- matrix(FALSE, p, p)
  for (j in 1:p) {
    r <- abs(R[j, ]); r[j] <- -Inf
    cand <- order(r, decreasing = TRUE)[1:dn]
    fit <- glmnet(X[, cand, drop = FALSE], X[, j], intercept = FALSE,
                  standardize = TRUE, nlambda = 50)
    rss <- colSums((X[, j] - predict(fit, X[, cand, drop = FALSE]))^2)
    ebic <- n * log(rss / n) + fit$df * (log(n) + 2 * log(dn))
    b <- as.numeric(coef(fit, s = fit$lambda[which.min(ebic)]))[-1]
    B[j, cand[b != 0]] <- TRUE
  }
  B & t(B)
}
## huge() fits the glasso path to cor(X), so the EBIC likelihood is evaluated
## at cor(X); the argument S is kept for compatibility and no longer used.
## Note: the data analysis uses lambda.min.ratio = 0.05 for the glasso and MB
## paths (0.01 in the simulation); the manuscript should say so.
glasso_ebic <- function(X, S, nlambda = 30, gamma = 0.5) {
  n <- nrow(X); p <- ncol(X)
  S <- cor(X)
  h <- huge(X, method = "glasso", nlambda = nlambda, lambda.min.ratio = 0.05,
            verbose = FALSE)
  eb <- sapply(seq_along(h$lambda), function(k) {
    Om <- as.matrix(h$icov[[k]])
    ne <- (sum(abs(Om) > 1e-10) - p) / 2
    ld <- as.numeric(determinant(Om, logarithm = TRUE)$modulus)
    n * (sum(S * Om) - ld) + ne * log(n) + 4 * gamma * ne * log(p)
  })
  list(sel = as.matrix(h$path[[which.min(eb)]]) != 0,
       path = lapply(h$path, function(A) as.matrix(A) != 0))
}
## top-K edge set of a score matrix, and the member of a path of graphs whose
## size is closest to K
topK <- function(R, K) {
  s <- abs(R[upper.tri(R)]); A <- matrix(FALSE, nrow(R), ncol(R))
  A[upper.tri(A)] <- rank(-s, ties.method = "first") <= K
  A | t(A)
}
closestK <- function(path, K) {
  ne <- sapply(path, function(A) sum(A[upper.tri(A)]))
  path[[which.min(abs(ne - K))]]
}

## ---------------------------------------------------------------------------
## 2. Data, profiling and cleaning
## ---------------------------------------------------------------------------
data(stockdata, package = "huge")
P <- stockdata$data                       # 1258 x 452 closing prices
sector <- stockdata$info[, 2]
ticker <- stockdata$info[, 1]
stopifnot(!anyNA(P), all(P > 0), !anyDuplicated(ticker))
R0 <- diff(log(P))                        # 1257 x 452 log returns

## Closing prices are not adjusted for stock splits: log returns close to
## log(1/2), log(1/3), log(2/3), ... mark corporate actions, not price moves.
## Returns with |r| > 0.25 are treated as such events and set to zero (the
## median daily return is 0.0004); the count is recorded.
flag <- abs(R0) > 0.25
Rc <- R0; Rc[flag] <- 0
profile <- data.frame(
  quantity = c("trading days (returns)", "stocks", "sectors",
               "cells flagged as corporate actions",
               "share of cells flagged", "stocks with a flagged cell",
               "median excess kurtosis before cleaning",
               "median excess kurtosis after cleaning",
               "share of variance, first principal component"),
  value = c(nrow(R0), ncol(R0), length(unique(sector)), sum(flag),
            mean(flag), sum(colSums(flag) > 0),
            median(apply(R0, 2, function(x) mean((x - mean(x))^4) / var(x)^2 - 3)),
            median(apply(Rc, 2, function(x) mean((x - mean(x))^4) / var(x)^2 - 3)),
            { ev <- eigen(cor(Rc), TRUE, only.values = TRUE)$values
              ev[1] / sum(ev) }))
write.csv(profile, file.path(OUT, "results", "data_profile.csv"),
          row.names = FALSE)

## Within each analysis window every stock is mapped to normal scores
## Phi^{-1}(rank / (m + 1)) (a nonparanormal transformation, Liu et al., 2009),
## which removes the heavy tails without altering the rank dependence, and
## then column-centred.
normal_scores <- function(Y) {
  Z <- apply(Y, 2, function(x) qnorm(rank(x, ties.method = "average") /
                                        (length(x) + 1)))
  scale(Z, center = TRUE, scale = FALSE)
}

## five consecutive non-overlapping windows of 251 trading days (about one
## year each, p / n = 1.8) and the full sample (p / n = 0.36)
nw <- 251
windows <- c(lapply(1:5, function(w) ((w - 1) * nw + 1):(w * nw)),
             list(1:nrow(Rc)))
wname <- c(paste0("W", 1:5), "Full")

## sector baseline: share of all pairs that lie within a sector
same_sector <- outer(sector, sector, "==")
ut <- upper.tri(same_sector)
baseline <- mean(same_sector[ut])

## ---------------------------------------------------------------------------
## 3. Window-level analysis
## ---------------------------------------------------------------------------
grid <- exp(seq(log(1e-2), log(20), length.out = 30))
lam_path <- exp(seq(log(1e-3), log(100), length.out = 60))
kap_rows <- list(); sel_all <- list(); marg_rows <- list(); curve_rows <- list()
KS <- c(100, 500); matched <- list()

for (w in seq_along(windows)) {
  X <- normal_scores(Rc[windows[[w]], ])
  n <- nrow(X) - 1                          # effective size after centring
  p <- ncol(X)
  S <- crossprod(X) / n
  sp <- spec_S(S)

  lam <- cv_lambda(X, grid, type = "ridge")
  R_plain <- pcor_of(ridge_Q(sp, lam))
  lam_a <- cv_lambda(X, grid, type = "alt")
  R_alt <- pcor_of(altridge_Q(sp, lam_a))
  R_rich <- pcor_of(richardson_Q(sp, n, lam, 2))
  set.seed(1)                               # fixed seed per window for MPS
  mps <- mps_select(X)

  kap_rows[[w]] <- data.frame(
    window = wname[w], n = n, p = p, ratio = p / n, lambda_cv = lam,
    kappa_hat = kappa_hat(sp, lam, n), s_hat = 1 - r_hat(sp, lam, n),
    kappa_floor_hat = kappa_floor_hat(sp, n),
    n_edges_mps = sum(mps$A[upper.tri(mps$A)]))

  ## effective penalty along a penalty path
  curve_rows[[w]] <- data.frame(window = wname[w], lambda = lam_path,
    kappa = sapply(lam_path, kappa_hat, sp = sp, n = n))

  ## agreement of ridge partial correlations with marginal correlations
  C <- cov2cor(S)
  marg_rows[[w]] <- data.frame(window = wname[w], lambda = lam_path[seq(1, 60, 3)],
    kappa = sapply(lam_path[seq(1, 60, 3)], kappa_hat, sp = sp, n = n),
    spearman = sapply(lam_path[seq(1, 60, 3)], function(l)
      cor(abs(pcor_of(ridge_Q(sp, l))[ut]), abs(C[ut]), method = "spearman")))

  gl <- glasso_ebic(X, S)
  set.seed(2)
  hmb <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.05,
              verbose = FALSE)
  ## matched-size comparison: every method reports K edges
  matched[[wname[w]]] <- lapply(KS, function(K) list(
    RT = topK(R_plain, K), vWP = topK(R_alt, K), RE = topK(R_rich, K),
    GLASSO = closestK(gl$path, K),
    MB = closestK(lapply(hmb$path, function(A) as.matrix(A) != 0), K)))
  sel <- list(
    vWP_lFDR = lfdr_select(R_alt),
    RT_lFDR  = lfdr_select(R_plain),
    RE_lFDR  = lfdr_select(R_rich),
    RS_lasso = screen_refit(X, R_plain, dn = ceiling(nrow(X) / log(nrow(X)))),
    RT_MPS   = mps$A,
    GLASSO   = gl$sel,
    MB       = as.matrix(huge.select(hmb, criterion = "ric",
                                     verbose = FALSE)$refit) != 0)
  for (k in names(sel)) diag(sel[[k]]) <- FALSE
  sel_all[[wname[w]]] <- sel
  cat(sprintf("window %s done: lambda=%.3f kappa=%.3f floor=%.3f\n", wname[w],
              lam, kap_rows[[w]]$kappa_hat, kap_rows[[w]]$kappa_floor_hat))
}
kap_tab <- do.call(rbind, kap_rows)

## leakage factor q_kappa = kappa (theta_max - theta_min) / (1 + kappa theta_max)
## of Proposition 3.4(c), with the extreme eigenvalues of Sigma replaced by
## those of the full-sample covariance matrix (p / n = 0.36), evaluated at the
## effective penalty of the full-sample fit
Xf <- normal_scores(Rc)
sf <- eigen(crossprod(Xf) / (nrow(Xf) - 1), TRUE, only.values = TRUE)$values
kf <- kap_tab$kappa_hat[kap_tab$window == "Full"]
qf <- kf * (1 / min(sf) - 1 / max(sf)) / (1 + kf / min(sf))
write.csv(data.frame(sigma_max_hat = max(sf), sigma_min_hat = min(sf),
                     kappa_full = kf, q_hat = qf),
          file.path(OUT, "results", "data_leakage_factor.csv"), row.names = FALSE)
write.csv(kap_tab, file.path(OUT, "results", "data_kappa.csv"), row.names = FALSE)
write.csv(do.call(rbind, curve_rows), file.path(OUT, "results", "data_kappa_curve.csv"),
          row.names = FALSE)
marg_tab <- do.call(rbind, marg_rows)
write.csv(marg_tab, file.path(OUT, "results", "data_marginal.csv"), row.names = FALSE)
saveRDS(sel_all, file.path(OUT, "results", "data_selections.rds"))

## ---------------------------------------------------------------------------
## 4. Method comparison: size, sector concordance, reproducibility
## ---------------------------------------------------------------------------
jaccard <- function(A, B) {
  a <- A[ut]; b <- B[ut]
  if (sum(a | b) == 0) NA else sum(a & b) / sum(a | b)
}
methods <- names(sel_all[[1]])
comp <- do.call(rbind, lapply(methods, function(m) {
  ne <- sapply(wname[1:5], function(w) sum(sel_all[[w]][[m]][ut]))
  within <- sapply(wname[1:5], function(w) {
    a <- sel_all[[w]][[m]][ut]; if (sum(a) == 0) NA else mean(same_sector[ut][a]) })
  jac <- sapply(1:4, function(w) jaccard(sel_all[[wname[w]]][[m]],
                                         sel_all[[wname[w + 1]]][[m]]))
  full <- sel_all[["Full"]][[m]]
  data.frame(method = m, edges_mean = mean(ne), edges_min = min(ne),
             edges_max = max(ne), within_sector = mean(within, na.rm = TRUE),
             jaccard_consecutive = mean(jac, na.rm = TRUE),
             edges_full = sum(full[ut]),
             within_full = if (sum(full[ut]) == 0) NA else
               mean(same_sector[ut][full[ut]]))
}))
comp$baseline <- baseline

## matched-size comparison (averages over the five yearly windows)
mt <- do.call(rbind, lapply(seq_along(KS), function(k) {
  do.call(rbind, lapply(names(matched[[1]][[k]]), function(m) {
    within <- sapply(wname[1:5], function(w) {
      a <- matched[[w]][[k]][[m]][ut]; mean(same_sector[ut][a]) })
    jac <- sapply(1:4, function(w) jaccard(matched[[wname[w]]][[k]][[m]],
                                           matched[[wname[w + 1]]][[k]][[m]]))
    size <- sapply(wname[1:5], function(w) sum(matched[[w]][[k]][[m]][ut]))
    data.frame(K = KS[k], method = m, mean_size = mean(size),
               within_sector = mean(within), jaccard_consecutive = mean(jac))
  }))
}))
write.csv(mt, file.path(OUT, "results", "data_matched.csv"), row.names = FALSE)
write.csv(comp, file.path(OUT, "results", "data_methods.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. Figures
## ---------------------------------------------------------------------------
curves <- do.call(rbind, curve_rows)
pdf(file.path(OUT, "figs", "fig_data_kappa.pdf"), width = 6.5, height = 3.4)
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1), mgp = c(2.4, 0.8, 0), cex = 0.85)
cols <- c(gray(seq(0.2, 0.65, length.out = 5)), "black")
plot(NA, xlim = range(lam_path), ylim = c(1e-3, 200), log = "xy",
     xlab = expression("nominal penalty " * lambda),
     ylab = expression("effective penalty " * hat(kappa)[lambda]))
for (w in seq_along(wname)) {
  cw <- curves[curves$window == wname[w], ]
  lines(cw$lambda, cw$kappa, col = cols[w], lty = if (w == 6) 2 else 1, lwd = 1.4)
}
abline(0, 1, col = "gray80")
legend("topleft", c("yearly windows (p/n = 1.8)", "full sample (p/n = 0.36)",
                    "identity"), lty = c(1, 2, 1),
       col = c(cols[3], "black", "gray80"), bty = "n", cex = 0.8)
plot(NA, xlim = range(marg_tab$kappa), ylim = c(0, 1), log = "x",
     xlab = expression("effective penalty " * hat(kappa)[lambda]),
     ylab = "Spearman correlation")
for (w in seq_along(wname)) {
  mw <- marg_tab[marg_tab$window == wname[w], ]
  lines(mw$kappa, mw$spearman, col = cols[w], lty = if (w == 6) 2 else 1, lwd = 1.4)
}
dev.off()

pdf(file.path(OUT, "figs", "fig_data_sector.pdf"), width = 6.5, height = 3.2)
par(mfrow = c(1, 2), mar = c(6.5, 4, 1, 1), mgp = c(2.4, 0.8, 0), cex = 0.8)
lab <- c("vWP-lFDR", "RT-lFDR", "RE-lFDR", "RS", "RT-MPS", "GLASSO", "MB")
barplot(comp$within_sector, names.arg = lab, las = 2, ylim = c(0, 1),
        ylab = "within-sector share of edges", col = "gray70", border = NA)
abline(h = baseline, lty = 2)
barplot(comp$jaccard_consecutive, names.arg = lab, las = 2, ylim = c(0, 0.6),
        ylab = "Jaccard index, consecutive windows", col = "gray70", border = NA)
dev.off()

## network of the calibrated ridge selection in the first window, coloured
## by sector; isolated stocks are omitted
A1 <- sel_all[["W1"]][["RT_MPS"]] * 1
keep <- which(rowSums(A1) > 0)
g1 <- graph_from_adjacency_matrix(A1[keep, keep], mode = "undirected")
sec_levels <- sort(unique(sector))
pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
         "#a6761d", "#666666", "#1f78b4", "#b2df8a")
vcol <- pal[match(sector[keep], sec_levels)]
set.seed(1)
pdf(file.path(OUT, "figs", "fig_data_network.pdf"), width = 6.5, height = 5)
par(mar = c(0, 0, 0, 0))
plot(g1, vertex.size = 3, vertex.label = NA, vertex.color = vcol,
     vertex.frame.color = NA, edge.color = "gray60", edge.width = 0.6,
     layout = layout_with_fr(g1))
legend("bottomleft", legend = sec_levels, pch = 19, col = pal, bty = "n",
       cex = 0.7)
dev.off()
## sector composition of the connected components of that network
cm <- components(g1)
comp_tab <- do.call(rbind, lapply(seq_len(cm$no), function(k) {
  v <- keep[cm$membership == k]
  data.frame(component = k, size = length(v),
             sectors = length(unique(sector[v])),
             main_sector = names(which.max(table(sector[v]))))
}))
comp_tab <- comp_tab[order(-comp_tab$size), ]
write.csv(comp_tab, file.path(OUT, "results", "data_components_W1.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()),
           file.path(OUT, "results", "sessionInfo_data.txt"))
cat("Data analysis finished.\n")
