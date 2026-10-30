###############################################################################
## Extensions of the real-data analysis for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Data   : `stockdata` from the R package huge, prepared exactly as in
##          data_stock_network.R (log returns, normal scores per window).
## Purpose: (1) the multi-penalty permutation selection (RT-MPS) with
##          B = 999 permutations, and its dependence on the random seed
##          against B = 99; (2) a factor-adjusted re-analysis (first 1, 2, 3
##          principal components removed, or the equal-weighted market return
##          regressed out); (3) block permutations as a check of the
##          exchangeability that MPS needs; (4) the Chebyshev leakage rate next
##          to q_kappa, and the downward bias of the smallest sample
##          eigenvalue; (5) why the Spearman correlation between ridge partial
##          correlations and sample correlations is close to zero at the
##          cross-validated penalty; (6) sensitivity to the cleaning of large
##          returns.
## Requirements: as data_stock_network.R, plus the base package parallel.
##          The building blocks (ridge estimators, cross-validation, MPS,
##          lFDR, screening, glasso with EBIC) are read from
##          data_stock_network.R, so both scripts run the same code.
## Usage  : Rscript data_extensions.R   (run from the code/ directory)
##          Each unit of work is cached in ../results/ext_cache/, so an
##          interrupted run resumes where it stopped.
## Output : ../results/data_mps_B999.csv, data_mps_B999_tables.csv,
##          data_factor.csv, data_factor_kappa.csv, data_blockperm.csv,
##          data_leakage_rates.csv, data_spearman.csv, data_cleaning.csv,
##          data_cleaning_sensitivity.csv, sessionInfo_extensions.txt
###############################################################################

if (nzchar(Sys.getenv("RIDGE_RLIB")))
  .libPaths(c(Sys.getenv("RIDGE_RLIB"), .libPaths()))
suppressPackageStartupMessages({
  library(huge)
  library(glmnet)
  library(fdrtool)
  library(parallel)
})
OUT <- normalizePath("..")
CACHE <- file.path(OUT, "results", "ext_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
NCORES <- max(1, min(23, detectCores() - 1))

## ---------------------------------------------------------------------------
## 1. Building blocks from the main data script (functions and constants only)
## ---------------------------------------------------------------------------
ex <- parse("data_stock_network.R")
is_def <- function(e) {
  if (!is.call(e) || !identical(e[[1]], as.name("<-")) || !is.name(e[[2]]))
    return(FALSE)
  (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) ||
    as.character(e[[2]]) %in% c("NPERM", "ALPHA", "MULT")
}
for (e in ex[vapply(ex, is_def, logical(1))]) eval(e, globalenv())
grid <- exp(seq(log(1e-2), log(20), length.out = 30))   # as in the main script
lam_path <- exp(seq(log(1e-3), log(100), length.out = 60))
KS <- c(100, 250, 500)

## ---------------------------------------------------------------------------
## 2. Data and cleaning rules
## ---------------------------------------------------------------------------
data(stockdata, package = "huge")
P <- stockdata$data
sector <- stockdata$info[, 2]
ticker <- stockdata$info[, 1]
R0 <- diff(log(P))
p <- ncol(R0)
big <- abs(R0) > 0.25
## log returns that match a split (1/2, 1/3, 2/3, 3/4, 1/4, 2/5) or a reverse
## split (5, 6, 8, 10) to within 0.05 are treated as corporate actions
ratios <- log(c(1/2, 1/3, 2/3, 3/4, 1/4, 2/5, 5, 6, 8, 10))
dist_ratio <- array(apply(abs(outer(R0, ratios, "-")), c(1, 2), min), dim(R0))
split_like <- big & dist_ratio < 0.05

## rule "zero"  : every |r| > 0.25 set to zero (main analysis)
## rule "winsor": every |r| > 0.25 kept, capped at +-0.25
## rule "split" : split-like returns set to zero, the others capped at +-0.25
clean_returns <- function(rule) {
  R <- R0
  if (rule == "zero") R[big] <- 0
  if (rule == "split") R[split_like] <- 0
  if (rule %in% c("winsor", "split")) R <- pmin(pmax(R, -0.25), 0.25)
  R
}
nw <- 251
windows <- c(lapply(1:5, function(w) ((w - 1) * nw + 1):(w * nw)),
             list(1:nrow(R0)))
wname <- c(paste0("W", 1:5), "Full")
same_sector <- outer(sector, sector, "==")
ut <- upper.tri(same_sector)
ss <- same_sector[ut]

## factor adjustment within a window. X: centred normal scores; Rw: cleaned
## returns of the window. Removing k principal components leaves n - 1 - k
## degrees of freedom; the residuals are rescaled to unit variance, so that
## ridge acts on the residual correlation matrix.
adjust_X <- function(X, Rw, how) {
  if (how == "none") return(list(X = X, n = nrow(X) - 1, k = 0))
  if (how == "mkt") {
    m <- rowMeans(Rw); m <- m - mean(m)
    Xr <- X - m %*% (crossprod(m, X) / sum(m^2)); k <- 1
  } else {
    k <- as.integer(sub("pc", "", how))
    s <- svd(X, nu = k, nv = k)
    Xr <- X - s$u %*% (s$d[1:k] * t(s$v))
  }
  n <- nrow(X) - 1 - k
  list(X = sweep(Xr, 2, sqrt(colSums(Xr^2) / n), "/"), n = n, k = k)
}
make_X <- function(rule, w, how = "none") {
  Rc <- clean_returns(rule)
  adjust_X(normal_scores(Rc[windows[[w]], ]), Rc[windows[[w]], ], how)
}

## ---------------------------------------------------------------------------
## 3. Units of work
## ---------------------------------------------------------------------------
## permutation of every column in blocks of L consecutive days, after a random
## circular shift (L = 1 is the ordinary permutation of mps_select)
block_perm <- function(X, L) {
  n <- nrow(X); b <- split(seq_len(n), ceiling(seq_len(n) / L))
  apply(X, 2, function(x) {
    x <- x[(seq_len(n) + sample.int(n, 1) - 2) %% n + 1]
    x[unlist(b[sample.int(length(b))], use.names = FALSE)]
  })
}
## the computation of mps_select, returning the permutation maxima as well;
## with L = 1 it uses the random numbers exactly as mps_select does
mps_run <- function(X, B, seed, L = 1, alpha = ALPHA, mult = MULT) {
  set.seed(seed)
  n <- nrow(X) - 1; p <- ncol(X); ut <- upper.tri(diag(p))
  S <- crossprod(X) / n
  lams <- mult * sum(diag(S)) / p
  sp <- spec_S(S)
  Robs <- lapply(lams, function(l) pcor_of(ridge_Q(sp, l)))
  Tobs <- sapply(Robs, function(R) max(abs(R[ut])))
  Tb <- matrix(0, B, length(lams))
  for (b in 1:B) {
    Xp <- if (L == 1) apply(X, 2, sample) else block_perm(X, L)
    Xb <- scale(Xp, center = TRUE, scale = FALSE)
    spb <- spec_S(crossprod(Xb) / n)
    Tb[b, ] <- sapply(lams, function(l) max(abs(pcor_of(ridge_Q(spb, l))[ut])))
  }
  k <- ceiling((1 - alpha / length(lams)) * (B + 1))
  tau <- sapply(seq_along(lams), function(g) sort(c(Tobs[g], Tb[, g]))[k])
  A <- matrix(FALSE, p, p)
  for (g in seq_along(lams)) A <- A | (abs(Robs[[g]]) > tau[g])
  list(Tobs = Tobs, Tb = Tb, tau = tau, lambda = lams, edges = which(A[ut]))
}
## everything else of the window analysis: penalties, effective penalty,
## rankings at matched sizes and (optionally) the selections of Table 8
fit_window <- function(X, n, seed, paths = TRUE, select = FALSE) {
  set.seed(seed)
  ut <- upper.tri(diag(ncol(X)))
  S <- crossprod(X) / n; sp <- spec_S(S)
  lam <- cv_lambda(X, grid, type = "ridge")
  lam_a <- cv_lambda(X, grid, type = "alt")
  R_plain <- pcor_of(ridge_Q(sp, lam))
  R_alt <- pcor_of(altridge_Q(sp, lam_a))
  R_rich <- pcor_of(richardson_Q(sp, n, lam, 2))
  C <- cov2cor(S)
  top <- function(R, K) which(topK(R, K)[ut])
  res <- list(lambda_cv = lam, lambda_alt = lam_a,
              kappa_hat = kappa_hat(sp, lam, n), s_hat = 1 - r_hat(sp, lam, n),
              kappa_floor_hat = kappa_floor_hat(sp, n),
              sigma_max = max(sp$val), spearman_cv =
                cor(abs(R_plain[ut]), abs(C[ut]), method = "spearman"),
              matched = lapply(KS, function(K) list(
                RT = top(R_plain, K), vWP = top(R_alt, K), RE = top(R_rich, K),
                COR = top(C, K))))
  if (paths) {
    gl <- glasso_ebic(X, S)
    hmb <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.05,
                verbose = FALSE)
    mbp <- lapply(hmb$path, function(A) as.matrix(A) != 0)
    for (k in seq_along(KS)) {
      res$matched[[k]]$GLASSO <- which(closestK(gl$path, KS[k])[ut])
      res$matched[[k]]$MB <- which(closestK(mbp, KS[k])[ut])
    }
    if (select) {
      safe <- function(R) tryCatch(which(lfdr_select(R)[ut]),
                                   error = function(e) integer(0))
      res$sel <- list(
        vWP_lFDR = safe(R_alt), RT_lFDR = safe(R_plain), RE_lFDR = safe(R_rich),
        RS_lasso = which(screen_refit(X, R_plain,
                           dn = ceiling(nrow(X) / log(nrow(X))))[ut]),
        GLASSO = which(gl$sel[ut]),
        MB = which((as.matrix(huge.select(hmb, criterion = "ric",
                                          verbose = FALSE)$refit) != 0)[ut]))
    }
  }
  res
}
run_job <- function(job) {
  f <- file.path(CACHE, paste0(job$id, ".rds"))
  if (file.exists(f)) return(job$id)
  t0 <- proc.time()[3]
  d <- make_X(job$rule, job$w, job$how)
  res <- if (job$type == "mps") mps_run(d$X, job$B, job$seed, job$L)
         else fit_window(d$X, d$n, 20260930 + job$w,
                         paths = job$w <= 5 || job$how == "none",
                         select = job$how == "none")
  res$seconds <- proc.time()[3] - t0
  saveRDS(res, f)
  job$id
}

jobs <- list()
add <- function(...) { j <- list(...); jobs[[j$id]] <<- j }
for (w in c(6, 1:5)) {
  for (s in 1:5) add(id = sprintf("mps_zero_%s_B999_s%d_L1", wname[w], s),
                     type = "mps", rule = "zero", w = w, how = "none",
                     B = 999, seed = s, L = 1)
  for (L in c(5, 20, 60)) add(id = sprintf("mps_zero_%s_B999_s1_L%d", wname[w], L),
                              type = "mps", rule = "zero", w = w, how = "none",
                              B = 999, seed = 1, L = L)
  for (r in c("winsor", "split")) add(id = sprintf("mps_%s_%s_B999_s1_L1", r, wname[w]),
                                      type = "mps", rule = r, w = w, how = "none",
                                      B = 999, seed = 1, L = 1)
}
for (w in c(6, 1:5)) {
  for (r in c("zero", "winsor", "split"))
    add(id = sprintf("fit_%s_%s_none", r, wname[w]), type = "fit", rule = r,
        w = w, how = "none")
  for (h in c("pc1", "pc2", "pc3", "mkt"))
    add(id = sprintf("fit_zero_%s_%s", wname[w], h), type = "fit", rule = "zero",
        w = w, how = h)
  for (s in 1:5) add(id = sprintf("mps_zero_%s_B99_s%d_L1", wname[w], s),
                     type = "mps", rule = "zero", w = w, how = "none",
                     B = 99, seed = s, L = 1)
}
todo <- Filter(function(j) !file.exists(file.path(CACHE, paste0(j$id, ".rds"))),
               jobs)
cat(length(jobs), "units of work,", length(todo), "to run on", NCORES, "cores\n")
t_start <- proc.time()[3]
if (length(todo)) {
  cl <- makePSOCKcluster(min(NCORES, length(todo)))
  clusterCall(cl, function(lp) { .libPaths(lp)
    suppressPackageStartupMessages({ library(huge); library(glmnet)
      library(fdrtool) }); NULL }, .libPaths())
  clusterExport(cl, setdiff(ls(globalenv()), c("cl", "todo", "jobs", "ex")))
  invisible(clusterApplyLB(cl, todo, run_job))
  stopCluster(cl)
}
cat(sprintf("units of work finished in %.1f minutes\n",
            (proc.time()[3] - t_start) / 60))
get <- function(id) readRDS(file.path(CACHE, paste0(id, ".rds")))

## ---------------------------------------------------------------------------
## 4. Task 1: RT-MPS with B = 999 against B = 99
## ---------------------------------------------------------------------------
within_of <- function(idx) if (length(idx)) mean(ss[idx]) else NA
jac <- function(a, b) { u <- length(union(a, b)); if (u) length(intersect(a, b)) / u else NA }
mps_rows <- list()
for (w in 1:6) for (B in c(99, 999)) for (s in 1:5) {
  m <- get(sprintf("mps_zero_%s_B%d_s%d_L1", wname[w], B, s))
  mps_rows[[length(mps_rows) + 1]] <- data.frame(window = wname[w], B = B,
    seed = as.character(s), edges = length(m$edges),
    within_sector = within_of(m$edges), tau1 = m$tau[1], tau2 = m$tau[2],
    tau3 = m$tau[3], tau4 = m$tau[4], minutes = m$seconds / 60)
}
## pooled threshold: the 5 x 999 permutation maxima of one window together
mps_pooled <- list()
for (w in 1:6) {
  ms <- lapply(1:5, function(s) get(sprintf("mps_zero_%s_B999_s%d_L1", wname[w], s)))
  Tb <- do.call(rbind, lapply(ms, `[[`, "Tb")); Tobs <- ms[[1]]$Tobs
  k <- ceiling((1 - ALPHA / length(MULT)) * (nrow(Tb) + 1))
  tau <- sapply(seq_along(MULT), function(g) sort(c(Tobs[g], Tb[, g]))[k])
  X <- make_X("zero", w)$X; n <- nrow(X) - 1
  sp <- spec_S(crossprod(X) / n)
  A <- matrix(FALSE, p, p)
  for (g in seq_along(MULT)) A <- A | (abs(pcor_of(ridge_Q(sp, ms[[1]]$lambda[g]))) > tau[g])
  mps_pooled[[wname[w]]] <- which(A[ut])
  mps_rows[[length(mps_rows) + 1]] <- data.frame(window = wname[w], B = nrow(Tb),
    seed = "pooled 1-5", edges = length(mps_pooled[[w]]),
    within_sector = within_of(mps_pooled[[w]]), tau1 = tau[1], tau2 = tau[2],
    tau3 = tau[3], tau4 = tau[4], minutes = NA)
}
mps_tab <- do.call(rbind, mps_rows)
write.csv(mps_tab, file.path(OUT, "results", "data_mps_B999.csv"), row.names = FALSE)

## the RT-MPS entries of Tables 7 and 8 for each run
mps_summary <- function(edge_list) {
  ne <- sapply(edge_list[1:5], length)
  data.frame(W1 = ne[1], W2 = ne[2], W3 = ne[3], W4 = ne[4], W5 = ne[5],
             edges_mean = mean(ne),
             within_sector = mean(sapply(edge_list[1:5], within_of), na.rm = TRUE),
             jaccard_consecutive = mean(sapply(1:4, function(w)
               jac(edge_list[[w]], edge_list[[w + 1]])), na.rm = TRUE),
             edges_full = length(edge_list[[6]]),
             within_full = within_of(edge_list[[6]]))
}
tab78 <- list()
for (B in c(99, 999)) for (s in 1:5) {
  el <- lapply(1:6, function(w) get(sprintf("mps_zero_%s_B%d_s%d_L1", wname[w], B, s))$edges)
  tab78[[length(tab78) + 1]] <- cbind(B = B, seed = as.character(s), mps_summary(el))
}
tab78[[length(tab78) + 1]] <- cbind(B = 4995, seed = "pooled 1-5",
                                    mps_summary(mps_pooled))
tab78 <- do.call(rbind, tab78)
write.csv(tab78, file.path(OUT, "results", "data_mps_B999_tables.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. Task 3: block permutations
## ---------------------------------------------------------------------------
acf1 <- function(x) cor(x[-1], x[-length(x)])
bp_rows <- list()
for (w in 1:6) {
  X <- make_X("zero", w)$X
  dep <- c(acf_level = median(apply(X, 2, acf1)),
           acf_square = median(apply(X^2, 2, acf1)))
  sd_tau <- apply(sapply(1:5, function(s)
    get(sprintf("mps_zero_%s_B999_s%d_L1", wname[w], s))$tau), 1, sd)
  for (L in c(1, 5, 20, 60)) {
    m <- get(sprintf("mps_zero_%s_B999_s1_L%d", wname[w], L))
    bp_rows[[length(bp_rows) + 1]] <- data.frame(window = wname[w], block = L,
      B = 999, seed = 1, edges = length(m$edges), within_sector = within_of(m$edges),
      tau1 = m$tau[1], tau2 = m$tau[2], tau3 = m$tau[3], tau4 = m$tau[4],
      null_mean_max1 = mean(m$Tb[, 1]), null_mean_max4 = mean(m$Tb[, 4]),
      sd_tau1_over_seeds = sd_tau[1], sd_tau4_over_seeds = sd_tau[4],
      median_acf1_level = dep[1], median_acf1_square = dep[2])
  }
}
write.csv(do.call(rbind, bp_rows), file.path(OUT, "results", "data_blockperm.csv"),
          row.names = FALSE)
## the RT-MPS entries of Tables 7 and 8 for each block length (seed 1)
bp_tab <- do.call(rbind, lapply(c(1, 5, 20, 60), function(L) cbind(block = L,
  mps_summary(lapply(1:6, function(w)
    get(sprintf("mps_zero_%s_B999_s1_L%d", wname[w], L))$edges)))))
write.csv(bp_tab, file.path(OUT, "results", "data_blockperm_tables.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 6. Task 2: factor-adjusted re-analysis
## ---------------------------------------------------------------------------
hows <- c("none", "pc1", "pc2", "pc3", "mkt")
fit_id <- function(w, h) if (h == "none") sprintf("fit_zero_%s_none", wname[w]) else
  sprintf("fit_zero_%s_%s", wname[w], h)
## leakage factor and Chebyshev rate with the extreme eigenvalues of a
## full-sample covariance matrix standing in for those of Sigma
q_rates <- function(kappa, smax, smin) {
  tmax <- 1 / smin; tmin <- 1 / smax
  r <- (1 + kappa * tmax) / (1 + kappa * tmin)
  c(q_kappa = kappa * (tmax - tmin) / (1 + kappa * tmax),
    q_C = (sqrt(r) - 1) / (sqrt(r) + 1),
    q_cond = { c0 <- (smax + kappa) / (smin + kappa); (sqrt(c0) - 1) / (sqrt(c0) + 1) })
}
full_spec <- lapply(hows, function(h) {
  d <- make_X("zero", 6, h)
  ev <- eigen(crossprod(d$X) / d$n, TRUE, only.values = TRUE)$values
  ev[seq_len(length(ev) - d$k)]          # drop the k zeros created by removing PCs
})
names(full_spec) <- hows
fk_rows <- list(); fm_rows <- list()
for (h in hows) {
  fits <- lapply(1:6, function(w) get(fit_id(w, h)))
  sv <- full_spec[[h]]
  for (w in 1:6) {
    f <- fits[[w]]
    fk_rows[[length(fk_rows) + 1]] <- data.frame(adjustment = h, window = wname[w],
      lambda_cv = f$lambda_cv, kappa_hat = f$kappa_hat, s_hat = f$s_hat,
      kappa_floor_hat = f$kappa_floor_hat, sigma_max_window = f$sigma_max,
      sigma_max_full = max(sv), sigma_min_full = min(sv),
      t(q_rates(f$kappa_hat, max(sv), min(sv))), spearman_cv = f$spearman_cv)
  }
  for (k in seq_along(KS)) for (m in names(fits[[1]]$matched[[k]])) {
    sets <- lapply(1:5, function(w) fits[[w]]$matched[[k]][[m]])
    fm_rows[[length(fm_rows) + 1]] <- data.frame(adjustment = h, K = KS[k],
      method = m, mean_size = mean(sapply(sets, length)),
      within_sector = mean(sapply(sets, within_of)),
      jaccard_consecutive = mean(sapply(1:4, function(w) jac(sets[[w]], sets[[w + 1]]))))
  }
}
write.csv(do.call(rbind, fm_rows), file.path(OUT, "results", "data_factor.csv"),
          row.names = FALSE)
fk_tab <- do.call(rbind, fk_rows)
write.csv(fk_tab, file.path(OUT, "results", "data_factor_kappa.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 7. Task 4: leakage rates and the bias of the extreme sample eigenvalues
## ---------------------------------------------------------------------------
## parametric bootstrap: if Sigma were the full-sample covariance matrix, how
## far would the extreme eigenvalues of a sample of the same size fall from it?
set.seed(4)
eig_bias <- t(sapply(hows, function(h) {
  d <- make_X("zero", 6, h); S <- crossprod(d$X) / d$n
  e <- eigen(S, TRUE); keep <- seq_len(p - d$k)
  root <- e$vectors[, keep] %*% (sqrt(e$values[keep]) * t(e$vectors[, keep]))
  bs <- replicate(20, {
    Z <- matrix(rnorm((d$n + 1) * p), d$n + 1, p) %*% root
    Z <- scale(Z, TRUE, FALSE)
    v <- eigen(crossprod(Z) / d$n, TRUE, only.values = TRUE)$values[keep]
    c(max(v), min(v))
  })
  c(ratio_max = mean(bs[1, ]) / max(e$values[keep]),
    ratio_min = mean(bs[2, ]) / min(e$values[keep]))
}))
lr_rows <- list()
for (h in hows) {
  sv <- full_spec[[h]]
  smin_c <- min(sv) / eig_bias[h, "ratio_min"]; smax_c <- max(sv) / eig_bias[h, "ratio_max"]
  for (w in 1:6) {
    k <- fk_tab$kappa_hat[fk_tab$adjustment == h & fk_tab$window == wname[w]]
    kf <- fk_tab$kappa_floor_hat[fk_tab$adjustment == h & fk_tab$window == wname[w]]
    a <- q_rates(k, max(sv), min(sv)); b <- q_rates(k, smax_c, smin_c)
    lr_rows[[length(lr_rows) + 1]] <- data.frame(adjustment = h, window = wname[w],
      kappa_hat = k, kappa_floor_hat = kf, sigma_max = max(sv), sigma_min = min(sv),
      q_kappa = a[1], q_C = a[2], q_cond_SigmaPlusKappa = a[3],
      ratio_q6 = (a[1] / a[2])^6,
      bias_ratio_min = eig_bias[h, "ratio_min"], bias_ratio_max = eig_bias[h, "ratio_max"],
      sigma_min_corrected = smin_c, sigma_max_corrected = smax_c,
      q_kappa_corrected = b[1], q_C_corrected = b[2])
  }
}
write.csv(do.call(rbind, lr_rows), file.path(OUT, "results", "data_leakage_rates.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## 8. Task 5: the Spearman correlation at the cross-validated penalty
## ---------------------------------------------------------------------------
sp_rows <- list()
set.seed(5)
for (w in 1:6) {
  X <- make_X("zero", w)$X; n <- nrow(X) - 1
  S <- crossprod(X) / n; sp <- spec_S(S); C <- cov2cor(S)
  lam <- get(fit_id(w, "none"))$lambda_cv
  rho <- abs(pcor_of(ridge_Q(sp, lam))[ut]); cabs <- abs(C[ut])
  e1 <- eigen(C, TRUE); l1 <- e1$vectors[, 1] * sqrt(e1$values[1])
  load <- abs(outer(l1, l1)[ut])
  res1 <- abs(cov2cor(crossprod(make_X("zero", w, "pc1")$X))[ut])
  res3 <- abs(cov2cor(crossprod(make_X("zero", w, "pc3")$X))[ut])
  Xn <- scale(apply(X, 2, sample), TRUE, FALSE)  # independent columns, same margins
  spn <- spec_S(crossprod(Xn) / n)
  rho_n <- abs(pcor_of(ridge_Q(spn, lam))[ut]); c_n <- abs(cov2cor(crossprod(Xn))[ut])
  rho_big <- abs(pcor_of(ridge_Q(sp, 100))[ut])
  sp_ <- function(a, b, i = TRUE) cor(a[i], b[i], method = "spearman")
  top <- rank(-rho) <= 1000
  sp_rows[[w]] <- data.frame(window = wname[w], lambda_cv = lam,
    kappa_hat = kappa_hat(sp, lam, n), kappa_over_sigma_max = kappa_hat(sp, lam, n) / max(sp$val),
    sp_rho_cor = sp_(rho, cabs),
    sp_cor_loadings = sp_(cabs, load),
    sp_rho_loadings = sp_(rho, load),
    sp_rho_resid1 = sp_(rho, res1),
    sp_rho_resid3 = sp_(rho, res3),
    sp_rho_cor_within = sp_(rho, cabs, ss),
    sp_rho_cor_between = sp_(rho, cabs, !ss),
    sp_rho_cor_null = sp_(rho_n, c_n),
    sp_rho_cor_lambda100 = sp_(rho_big, cabs),
    share_neg_rho = mean(pcor_of(ridge_Q(sp, lam))[ut] < 0),
    top1000_mean_cor_rank_pct = mean(rank(cabs)[top]) / length(cabs),
    top1000_within = mean(ss[top]),
    mean_abs_cor = mean(cabs), sd_abs_cor = sd(cabs), share_cor_positive = mean(C[ut] > 0))
}
write.csv(do.call(rbind, sp_rows), file.path(OUT, "results", "data_spearman.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## 9. Task 6: cleaning of large returns
## ---------------------------------------------------------------------------
win_of <- function(t) ifelse(t > 5 * nw, "after W5", paste0("W", pmin(ceiling(t / nw), 5)))
ix <- which(big, arr.ind = TRUE)
cl_tab <- as.data.frame.matrix(table(window = win_of(ix[, 1]),
  class = ifelse(split_like[ix], "split_like",
                 ifelse(R0[ix] > 0, "other_positive", "other_negative"))))
cl_tab <- cbind(window = rownames(cl_tab), cl_tab, total = rowSums(cl_tab))
write.csv(cl_tab, file.path(OUT, "results", "data_cleaning.csv"), row.names = FALSE)
write.csv(data.frame(row = ix[, 1], window = win_of(ix[, 1]), ticker = ticker[ix[, 2]],
                     sector = sector[ix[, 2]], return = R0[ix],
                     distance_to_split_ratio = dist_ratio[ix], split_like = split_like[ix]),
          file.path(OUT, "results", "data_cleaning_cells.csv"), row.names = FALSE)

sens <- list()
for (r in c("zero", "winsor", "split")) {
  fits <- lapply(1:6, function(w) get(sprintf("fit_%s_%s_none", r, wname[w])))
  mps <- lapply(1:6, function(w) get(if (r == "zero")
    sprintf("mps_zero_%s_B999_s1_L1", wname[w]) else
    sprintf("mps_%s_%s_B999_s1_L1", r, wname[w]))$edges)
  for (w in 1:6) sens[[length(sens) + 1]] <- data.frame(rule = r, table = "7",
    window = wname[w], method = "", lambda_cv = fits[[w]]$lambda_cv,
    kappa_hat = fits[[w]]$kappa_hat, s_hat = fits[[w]]$s_hat,
    kappa_floor_hat = fits[[w]]$kappa_floor_hat, edges = length(mps[[w]]),
    within_sector = NA, jaccard_consecutive = NA, edges_full = NA, within_full = NA)
  meth <- c(names(fits[[1]]$sel)[1:4], "RT_MPS", "GLASSO", "MB")
  for (m in meth) {
    el <- if (m == "RT_MPS") mps else lapply(fits, function(f) f$sel[[m]])
    s8 <- mps_summary(el)
    sens[[length(sens) + 1]] <- data.frame(rule = r, table = "8", window = "",
      method = m, lambda_cv = NA, kappa_hat = NA, s_hat = NA, kappa_floor_hat = NA,
      edges = s8$edges_mean, within_sector = s8$within_sector,
      jaccard_consecutive = s8$jaccard_consecutive, edges_full = s8$edges_full,
      within_full = s8$within_full)
  }
  for (k in c(1, 3)) for (m in c("RT", "vWP", "RE", "GLASSO", "MB")) {
    sets <- lapply(1:5, function(w) fits[[w]]$matched[[k]][[m]])
    sens[[length(sens) + 1]] <- data.frame(rule = r, table = paste0("9 K=", KS[k]),
      window = "", method = m, lambda_cv = NA, kappa_hat = NA, s_hat = NA,
      kappa_floor_hat = NA, edges = mean(sapply(sets, length)),
      within_sector = mean(sapply(sets, within_of)),
      jaccard_consecutive = mean(sapply(1:4, function(w) jac(sets[[w]], sets[[w + 1]]))),
      edges_full = NA, within_full = NA)
  }
}
write.csv(do.call(rbind, sens), file.path(OUT, "results", "data_cleaning_sensitivity.csv"),
          row.names = FALSE)

writeLines(c(sprintf("Total run time of the units of work: %.1f minutes on %d cores",
                     (proc.time()[3] - t_start) / 60, NCORES),
             capture.output(sessionInfo())),
           file.path(OUT, "results", "sessionInfo_extensions.txt"))
cat("Extensions finished.\n")
