###############################################################################
## Further analyses of the stock return data for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Data   : `stockdata` from the R package huge, prepared as in
##          data_stock_network.R (log returns, |r| > 0.25 set to zero, normal
##          scores within each window).
## Purpose: (1) networks compared at a matched nominal penalty and at a
##          matched estimated effective penalty: the full sample against each
##          yearly window, and random subsamples of several sizes; a
##          data-based leakage diagnostic for each window;
##          (2) the graphical lasso and neighbourhood selection with other
##          tuning rules (finer paths, StARS, EBIC with gamma = 0.5);
##          (3) reproducibility of every method at matched network sizes,
##          with GeneNet and beam added;
##          (4) the effective penalty and the network comparison after
##          GARCH(1,1) standardization of each return series and after
##          removal of the leading principal components, and a semi-synthetic
##          check with the full-sample correlation matrix as Sigma
##          (Gaussian rows, multivariate t rows, t rows with normal scores).
## Requirements: R >= 4.1; packages huge, glmnet, fdrtool, GeneNet, beam,
##          fGarch, parallel. The ridge building blocks are read from
##          data_stock_network.R and data_extensions.R, so all scripts run
##          the same code.
## Usage  : Rscript data_checks.R   (run from the code/ directory)
##          Each unit of work is saved in ../results/checks_cache/, so an
##          interrupted run resumes where it stopped.
## Output : ../results/data_checks_*.csv, sessionInfo_checks.txt
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
CACHE <- file.path(OUT, "results", "checks_cache")
EXT_CACHE <- file.path(OUT, "results", "ext_cache")
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
NCORES <- as.integer(Sys.getenv("RIDGE_CORES", "4"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

## ---------------------------------------------------------------------------
## 1. Building blocks from the two earlier data scripts
## ---------------------------------------------------------------------------
is_def <- function(e, keep = NULL) {
  if (!is.call(e) || !identical(e[[1]], as.name("<-")) || !is.name(e[[2]]))
    return(FALSE)
  nm <- as.character(e[[2]])
  if (!is.null(keep)) return(nm %in% keep)
  (is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) ||
    nm %in% c("NPERM", "ALPHA", "MULT")
}
ex <- parse("data_stock_network.R")
for (e in ex[vapply(ex, is_def, logical(1))]) eval(e, globalenv())
ex2 <- parse("data_extensions.R")
for (e in ex2[vapply(ex2, is_def, logical(1),
                     keep = c("adjust_X", "mps_run", "block_perm", "q_rates"))])
  eval(e, globalenv())
rm(ex, ex2, e)
grid <- exp(seq(log(1e-2), log(20), length.out = 30))   # as in the main script
NTOP <- 3000                     # length of the stored head of each ranking
KS <- c(100, 250, 500)

## ---------------------------------------------------------------------------
## 2. Data
## ---------------------------------------------------------------------------
data(stockdata, package = "huge")
P <- stockdata$data
sector <- stockdata$info[, 2]
R0 <- diff(log(P))
Rc <- R0; Rc[abs(R0) > 0.25] <- 0           # cleaning of the main analysis
p <- ncol(Rc); nT <- nrow(Rc)
nw <- 251
windows <- c(lapply(1:5, function(w) ((w - 1) * nw + 1):(w * nw)),
             list(1:nT))
wname <- c(paste0("W", 1:5), "Full")
same_sector <- outer(sector, sector, "==")
ut <- upper.tri(same_sector)
ss <- same_sector[ut]
m_pairs <- sum(ut)
pair_index <- matrix(0L, p, p); pair_index[ut] <- seq_len(m_pairs)

## indices (in upper-triangle order) of the |score| ranking, best first
rank_head <- function(R, len = NTOP) {
  s <- abs(R[ut]); order(-s)[seq_len(len)]
}
edges_of <- function(A) which(A[ut])

## ---------------------------------------------------------------------------
## 3. GARCH(1,1) standardization of each return series
## ---------------------------------------------------------------------------
## Each cleaned series (in per cent) is fitted by a GARCH(1,1) model with
## Gaussian innovations over the whole period; the standardized residuals
## replace the returns. If a fit fails, or is not stationary, the series is
## divided by an exponentially weighted volatility (decay 0.94) instead.
garch_one <- function(x) {
  x <- 100 * x
  fit <- tryCatch(suppressWarnings(fGarch::garchFit(~ garch(1, 1), data = x,
                  cond.dist = "norm", include.mean = TRUE, trace = FALSE)),
                  error = function(e) NULL)
  if (!is.null(fit)) {
    cf <- fit@fit$coef
    z <- as.numeric(fGarch::residuals(fit, standardize = TRUE))
    if (all(is.finite(z)) && isTRUE(cf["alpha1"] + cf["beta1"] < 1))
      return(list(z = z, coef = cf[c("mu", "omega", "alpha1", "beta1")],
                  method = "garch"))
  }
  v <- numeric(length(x)); v[1] <- var(x); xm <- x - mean(x)
  for (t in 2:length(x)) v[t] <- 0.94 * v[t - 1] + 0.06 * xm[t - 1]^2
  list(z = xm / sqrt(v), coef = c(mu = mean(x), omega = NA, alpha1 = NA,
                                  beta1 = NA), method = "ewma")
}

## design matrix of a window. version "orig": normal scores of the cleaned
## returns (the main analysis); "garch": normal scores of the GARCH residuals
Zg <- NULL
make_Xv <- function(ver, w, rows = windows[[w]]) {
  Y <- if (ver == "orig") Rc[rows, ] else Zg[rows, ]
  normal_scores(Y)
}

## ---------------------------------------------------------------------------
## 4. Units of work
## ---------------------------------------------------------------------------
## (a) RT-MPS with B = 999 and seed 1. For the original data this is the run
##     of data_extensions.R, which is reused when present.
job_mps <- function(ver, w) {
  f <- file.path(EXT_CACHE, sprintf("mps_zero_%s_B999_s1_L1.rds", wname[w]))
  if (ver == "orig" && file.exists(f)) return(readRDS(f))
  X <- make_Xv(ver, w)
  mps_run(X, 999, 1)
}

## (b) ridge fits, GeneNet, beam, selections, effective penalties and the
##     leakage diagnostic of one window
lfdr_safe <- function(R) tryCatch(edges_of(lfdr_select(R)),
                                  error = function(e) integer(0))
## ridge partial correlations at the nominal penalty whose estimated
## effective penalty is kappa (NULL if kappa is below the limit as lambda -> 0)
pcor_at_kappa <- function(sp, n, kappa) {
  if (kappa <= kappa_floor_hat(sp, n) * (1 + 1e-8)) return(NULL)
  l <- lambda_for_kappa(sp, n, kappa)
  list(lambda = l, R = pcor_of(ridge_Q(sp, l)))
}
## leakage diagnostic: change of the ridge partial correlations when the
## estimated effective penalty is halved (or doubled), on the same data
kappa_shift <- function(sp, n, R0_, kappa, factor) {
  a <- pcor_at_kappa(sp, n, factor * kappa)
  if (is.null(a)) return(c(lambda = NA, max = NA, q99 = NA, q999 = NA,
                           mean = NA, jac100 = NA, jac500 = NA, spearman = NA))
  d <- abs(a$R[ut] - R0_[ut])
  t0 <- rank_head(R0_, 500); t1 <- rank_head(a$R, 500)
  c(lambda = a$lambda, max = max(d), q99 = unname(quantile(d, 0.99)),
    q999 = unname(quantile(d, 0.999)), mean = mean(d),
    jac100 = length(intersect(t0[1:100], t1[1:100])) /
      length(union(t0[1:100], t1[1:100])),
    jac500 = length(intersect(t0, t1)) / length(union(t0, t1)),
    spearman = cor(abs(a$R[ut]), abs(R0_[ut]), method = "spearman"))
}
## estimated effective penalty when the k largest sample eigenvalues are
## left out at the same nominal penalty (k spikes removed, n - k degrees of
## freedom); one spike changes r_hat by at most 1/n
kappa_drop <- function(sp, lambda, n, k) {
  v <- sp$val[-seq_len(k)]
  lambda / (1 - sum(v / (v + lambda)) / (n - k))
}

job_core <- function(ver, w, mps) {
  set.seed(20260930 + w)              # the folds of data_extensions.R
  X <- make_Xv(ver, w); n <- nrow(X) - 1
  S <- crossprod(X) / n; sp <- spec_S(S)
  lam <- cv_lambda(X, grid, type = "ridge")
  lam_a <- cv_lambda(X, grid, type = "alt")
  k_hat <- kappa_hat(sp, lam, n)
  R_plain <- pcor_of(ridge_Q(sp, lam))
  R_alt <- pcor_of(altridge_Q(sp, lam_a))
  R_rich <- pcor_of(richardson_Q(sp, n, lam, 2))
  C <- cov2cor(S)
  ## RT-MPS score: largest ratio of |rho| to the permutation threshold
  mps_score <- Reduce(pmax, lapply(seq_along(mps$lambda), function(g)
    abs(pcor_of(ridge_Q(sp, mps$lambda[g]))) / mps$tau[g]))
  ## GeneNet (static shrinkage of the correlation matrix, analytic lambda)
  ## and beam (conjugate Bayes with an empirical-Bayes shrinkage alpha). Both
  ## are ridge inverses of the correlation matrix R with penalty lambda / (1 -
  ## lambda) and alpha / (1 - alpha); the estimated effective penalty is
  ## computed from the spectrum of R and reported in units of S (times the
  ## mean variance of the normal scores).
  Rm <- cor(X); spR <- spec_S(Rm); vbar <- mean(diag(S))
  gn <- GeneNet::ggm.estimate.pcor(X, method = "static", verbose = FALSE)
  l_gn <- attr(gn, "lambda"); gn <- unclass(gn); attributes(gn) <- list(dim = c(p, p))
  nt <- GeneNet::network.test.edges(gn, plot = FALSE, verbose = FALSE)
  gn_sel <- sort(pair_index[cbind(pmin(nt$node1, nt$node2),
                                  pmax(nt$node1, nt$node2))][nt$prob >= 0.8])
  bm <- beam::beam(X, type = "conditional", verbose = FALSE)
  a_bm <- bm@alphaOpt
  bm_cor <- matrix(0, p, p); bm_cor[ut] <- bm@table[, "p_cor"]; bm_cor <- bm_cor + t(bm_cor)
  bs <- beam::beam.select(bm, thres = 0.1, method = "BH")
  bm_sel <- sort(as.integer(rownames(bs@conditional)))
  ## checks that both are ridge inverses of R
  chk_gn <- max(abs(pcor_of(ridge_Q(spR, l_gn / (1 - l_gn))) - gn))
  chk_bm <- max(abs(pcor_of(ridge_Q(spR, a_bm / (1 - a_bm)))[ut] - bm_cor[ut]))
  ## factor-adjusted effective penalties: k = 1, 2, 3 principal components
  ## removed and residuals rescaled (as in data_extensions.R), with their own
  ## cross-validated penalty
  fac <- do.call(rbind, lapply(1:3, function(k) {
    d <- adjust_X(X, NULL, paste0("pc", k))
    set.seed(20260930 + w)
    lk <- cv_lambda(d$X, grid, type = "ridge")
    spk <- spec_S(crossprod(d$X) / d$n)
    spk$val <- spk$val[seq_len(p - k)]; spk$vec <- spk$vec[, seq_len(p - k)]
    data.frame(k = k, lambda_cv_adj = lk, kappa_hat_adj = kappa_hat(spk, lk, d$n),
               kappa_floor_adj = kappa_floor_hat(spk, d$n),
               kappa_adj_at_lambda = kappa_hat(spk, lam, d$n),
               kappa_drop_at_lambda = kappa_drop(sp, lam, n, k),
               sigma_max_adj = max(spk$val))
  }))
  list(
    n = n, lambda_cv = lam, lambda_alt = lam_a, kappa_hat = k_hat,
    s_hat = 1 - r_hat(sp, lam, n), kappa_floor_hat = kappa_floor_hat(sp, n),
    sigma_max = max(sp$val), trace_S = sum(sp$val), vbar = vbar,
    genenet_lambda = l_gn, genenet_ridge = l_gn / (1 - l_gn),
    genenet_kappa = vbar * kappa_hat(spR, l_gn / (1 - l_gn), n), genenet_check = chk_gn,
    beam_alpha = a_bm, beam_ridge = a_bm / (1 - a_bm),
    beam_kappa = vbar * kappa_hat(spR, a_bm / (1 - a_bm), n), beam_check = chk_bm,
    rt_kappa_corscale = vbar * kappa_hat(spR, lam / vbar, n),
    factor = fac,
    shift_half = kappa_shift(sp, n, R_plain, k_hat, 0.5),
    shift_double = kappa_shift(sp, n, R_plain, k_hat, 2),
    rank = list(RT = rank_head(R_plain), vWP = rank_head(R_alt),
                RE = rank_head(R_rich), MPS = rank_head(mps_score),
                GeneNet = rank_head(gn), beam = rank_head(bm_cor),
                COR = rank_head(C)),
    sel = list(vWP_lFDR = lfdr_safe(R_alt), RT_lFDR = lfdr_safe(R_plain),
               RE_lFDR = lfdr_safe(R_rich), RT_MPS = sort(mps$edges),
               GeneNet = gn_sel, beam = bm_sel))
}

## (c) the graphical lasso: EBIC (gamma = 0.5) on the path of the main
##     analysis and on a longer, finer path, and StARS
ebic_path <- function(h, X, gamma = 0.5) {
  n <- nrow(X); S <- cor(X)
  sapply(seq_along(h$lambda), function(k) {
    Om <- as.matrix(h$icov[[k]])
    ne <- (sum(abs(Om) > 1e-10) - p) / 2
    ld <- as.numeric(determinant(Om, logarithm = TRUE)$modulus)
    n * (sum(S * Om) - ld) + ne * log(n) + 4 * gamma * ne * log(p)
  })
}
path_edges <- function(path) lapply(path, function(A) edges_of(as.matrix(A) != 0))
job_glasso <- function(ver, w) {
  X <- make_Xv(ver, w)
  full <- w == 6
  res <- list()
  ## path of the main analysis (30 penalties down to 0.05 times the largest)
  ## in the full sample, where EBIC selected its last member; in the windows
  ## the selections of the main analysis are read from data_selections.rds
  if (full) {
    h0 <- huge(X, method = "glasso", nlambda = 30, lambda.min.ratio = 0.05,
               verbose = FALSE)
    e0 <- ebic_path(h0, X)
    res$old <- list(lambda = h0$lambda, ebic = e0, opt = which.min(e0),
                    sel = edges_of(as.matrix(h0$path[[which.min(e0)]]) != 0),
                    size = sapply(h0$path, function(A) sum(A) / 2))
    rm(h0); gc()
  }
  ## finer path of 60 penalties: down to 0.02 times the largest in the full
  ## sample, and to 0.1 in the windows, where EBIC selects near the top of
  ## the path and the matched sizes need at most a few thousand edges
  h1 <- huge(X, method = "glasso", nlambda = 60,
             lambda.min.ratio = if (full) 0.02 else 0.1, verbose = FALSE)
  e1 <- ebic_path(h1, X)
  res$fine <- list(lambda = h1$lambda, ebic = e1, opt = which.min(e1),
                   sel = edges_of(as.matrix(h1$path[[which.min(e1)]]) != 0),
                   path = path_edges(h1$path),
                   size = sapply(h1$path, function(A) sum(A) / 2))
  rm(h0, h1); gc()
  ## StARS (beta = 0.1, 20 subsamples of size 10 sqrt(n)) on 30 penalties
  ## down to 0.2 times the largest
  set.seed(7000 + w + if (ver == "garch") 100 else 0)
  h2 <- huge(X, method = "glasso", nlambda = 30, lambda.min.ratio = 0.2,
             verbose = FALSE)
  s2 <- huge.select(h2, criterion = "stars", stars.thresh = 0.1, rep.num = 20,
                    verbose = FALSE)
  res$stars <- list(lambda = h2$lambda, variability = s2$variability,
                    opt = s2$opt.index, sel = edges_of(as.matrix(s2$refit) != 0),
                    size = sapply(h2$path, function(A) sum(A) / 2))
  res
}

## (d) neighbourhood selection: RIC on the path of the main analysis and on
##     a finer path, StARS, and nodewise EBIC (gamma = 0.5)
mb_ebic <- function(X, gamma = 0.5) {
  n <- nrow(X); B <- matrix(FALSE, p, p)
  for (j in 1:p) {
    fit <- glmnet(X[, -j], X[, j], intercept = FALSE, standardize = TRUE,
                  nlambda = 100, lambda.min.ratio = 0.01)
    rss <- colSums((X[, j] - predict(fit, X[, -j]))^2)
    eb <- n * log(rss / n) + fit$df * (log(n) + 2 * gamma * log(p - 1))
    b <- as.numeric(coef(fit, s = fit$lambda[which.min(eb)]))[-1]
    B[j, (1:p)[-j][b != 0]] <- TRUE
  }
  list(or = edges_of(B | t(B)), and = edges_of(B & t(B)))
}
job_mb <- function(ver, w) {
  X <- make_Xv(ver, w)
  res <- list()
  set.seed(2)                         # as before the RIC call of the main script
  h0 <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.05,
             verbose = FALSE)
  r0 <- huge.select(h0, criterion = "ric", verbose = FALSE)
  res$old <- list(lambda = h0$lambda, opt = which.min(abs(h0$lambda - r0$opt.lambda)),
                  sel = edges_of(as.matrix(r0$refit) != 0),
                  size = sapply(h0$path, function(A) sum(A) / 2))
  h1 <- huge(X, method = "mb", nlambda = 60, lambda.min.ratio = 0.01,
             verbose = FALSE)
  set.seed(3000 + w + if (ver == "garch") 100 else 0)
  r1 <- huge.select(h1, criterion = "ric", verbose = FALSE)
  s1 <- huge.select(h1, criterion = "stars", stars.thresh = 0.1, rep.num = 20,
                    verbose = FALSE)
  res$ric_fine <- list(opt = which.min(abs(h1$lambda - r1$opt.lambda)), sel = edges_of(as.matrix(r1$refit) != 0))
  res$stars <- list(lambda = h1$lambda, variability = s1$variability,
                    opt = s1$opt.index, sel = edges_of(as.matrix(s1$refit) != 0))
  res$fine <- list(lambda = h1$lambda, path = path_edges(h1$path),
                   size = sapply(h1$path, function(A) sum(A) / 2))
  res$ebic <- mb_ebic(X)
  res
}

## (e) random subsamples of the whole period. Networks at the common nominal
##     penalty LAM0 and at the common estimated effective penalties KAP0.
LAM0 <- 0.392                    # cross-validated penalty of W1 and W5
KAP0 <- c(0.5, 1.0)
NSUB <- c(150, 251, 400, 600, 900, 1257)
job_sub <- function(nsub, r) {
  set.seed(100000 + 1000 * which(NSUB == nsub) + r)
  rows <- if (nsub == nT) 1:nT else sort(sample.int(nT, nsub))
  X <- make_Xv("orig", 6, rows); n <- nrow(X) - 1
  sp <- spec_S(crossprod(X) / n)
  one <- function(lambda) {
    if (is.na(lambda)) return(NULL)
    R <- pcor_of(ridge_Q(sp, lambda))
    list(lambda = lambda, kappa = kappa_hat(sp, lambda, n),
         lfdr = lfdr_safe(R), top = rank_head(R, 1000))
  }
  at <- c(list(lam = one(LAM0)), lapply(KAP0, function(k)
    one(if (k > kappa_floor_hat(sp, n)) lambda_for_kappa(sp, n, k) else NA)))
  names(at) <- c("lambda", paste0("kappa", KAP0))
  list(n = n, floor = kappa_floor_hat(sp, n), at = at)
}

## (f) semi-synthetic check: Sigma is the full-sample correlation matrix of
##     the normal scores (a spiked spectrum, first eigenvalue about 120),
##     n = 251 rows as in a yearly window
SIM_LAM <- c(0.05, 0.392, 1)
job_sim <- function(dist, chunk, reps = 10) {
  set.seed(500000 + 100 * match(dist, c("gauss", "t5", "t5ns")) + chunk)
  Xf <- make_Xv("orig", 6); Sig <- cov2cor(crossprod(Xf))
  e <- eigen(Sig, symmetric = TRUE)
  root <- e$vectors %*% (sqrt(pmax(e$values, 0)) * t(e$vectors))
  acc <- lapply(SIM_LAM, function(l) list(sum = numeric(m_pairs), sq = numeric(m_pairs)))
  kap <- matrix(NA, reps, length(SIM_LAM)); flo <- numeric(reps)
  shift <- matrix(NA, reps, 2)
  for (b in 1:reps) {
    Z <- matrix(rnorm(nw * p), nw, p)
    if (dist != "gauss") Z <- Z / sqrt(rchisq(nw, 5) / 3)  # t5 rows, unit variance
    X <- Z %*% root
    X <- if (dist == "t5ns") normal_scores(X) else scale(X, TRUE, FALSE)
    n <- nrow(X) - 1; sp <- spec_S(crossprod(X) / n)
    flo[b] <- kappa_floor_hat(sp, n)
    for (g in seq_along(SIM_LAM)) {
      R <- pcor_of(ridge_Q(sp, SIM_LAM[g]))[ut]
      acc[[g]]$sum <- acc[[g]]$sum + R; acc[[g]]$sq <- acc[[g]]$sq + R^2
      kap[b, g] <- kappa_hat(sp, SIM_LAM[g], n)
    }
    R1 <- pcor_of(ridge_Q(sp, 0.392))
    sh <- kappa_shift(sp, n, R1, kap[b, 2], 0.5)
    shift[b, ] <- sh[c("max", "q99")]
  }
  list(acc = acc, kappa = kap, floor = flo, shift = shift, reps = reps)
}

run_unit <- function(job) {
  f <- file.path(CACHE, paste0(job$id, ".rds"))
  if (file.exists(f)) return(job$id)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  t0 <- proc.time()[3]
  res <- tryCatch(switch(job$type,
    mps = job_mps(job$ver, job$w),
    core = job_core(job$ver, job$w,
                    readRDS(file.path(CACHE, sprintf("mps_%s_%s.rds", job$ver, wname[job$w])))),
    glasso = job_glasso(job$ver, job$w),
    mb = job_mb(job$ver, job$w),
    sub = job_sub(job$nsub, job$r),
    sim = job_sim(job$dist, job$chunk)),
    error = function(e) e)
  if (inherits(res, "error")) {
    cat(job$id, ":", conditionMessage(res), "\n",
        file = file.path(CACHE, "errors.log"), append = TRUE)
    return(paste(job$id, "failed"))
  }
  res$seconds <- proc.time()[3] - t0
  saveRDS(res, f)
  job$id
}

## ---------------------------------------------------------------------------
## 5. Running the units of work
## ---------------------------------------------------------------------------
t_start <- proc.time()[3]
gfile <- file.path(CACHE, "garch_residuals.rds")
if (!file.exists(gfile)) {
  gl <- lapply(seq_len(p), function(j) garch_one(Rc[, j]))   # about 5 minutes
  saveRDS(list(z = sapply(gl, `[[`, "z"), coef = t(sapply(gl, `[[`, "coef")),
               method = sapply(gl, `[[`, "method")), gfile)
}
gres <- readRDS(gfile)
Zg <- gres$z

jobs <- list()
add <- function(...) { j <- list(...); jobs[[j$id]] <<- j }
for (v in c("orig", "garch")) for (w in c(6, 1:5)) {
  add(id = sprintf("glasso_%s_%s", v, wname[w]), type = "glasso", ver = v, w = w)
  add(id = sprintf("mb_%s_%s", v, wname[w]), type = "mb", ver = v, w = w)
  add(id = sprintf("mps_%s_%s", v, wname[w]), type = "mps", ver = v, w = w)
}
for (d in c("gauss", "t5", "t5ns")) for (ch in 1:4)
  add(id = sprintf("sim_%s_%d", d, ch), type = "sim", dist = d, chunk = ch)
for (ns in NSUB) for (r in if (ns == nT) 1 else 1:20)
  add(id = sprintf("sub_%d_%02d", ns, r), type = "sub", nsub = ns, r = r)
core_jobs <- list()
for (v in c("orig", "garch")) for (w in c(6, 1:5))
  core_jobs[[length(core_jobs) + 1]] <- list(id = sprintf("core_%s_%s", v, wname[w]),
                                             type = "core", ver = v, w = w)
pending <- function(js) Filter(function(j)
  !file.exists(file.path(CACHE, paste0(j$id, ".rds"))), js)
todo <- c(pending(jobs), pending(core_jobs))
cat(length(jobs) + length(core_jobs), "units of work,", length(todo),
    "to run on", NCORES, "cores\n")
if (length(todo)) {
  cl <- makePSOCKcluster(min(NCORES, length(todo)))
  clusterCall(cl, function(lp) { .libPaths(lp)
    suppressPackageStartupMessages({ library(huge); library(glmnet)
      library(fdrtool) }); NULL }, .libPaths())
  clusterExport(cl, setdiff(ls(globalenv()), c("cl", "todo", "jobs", "core_jobs")))
  ## the MPS runs come first, because the core fits need their thresholds
  invisible(clusterApplyLB(cl, pending(jobs), run_unit))
  invisible(clusterApplyLB(cl, pending(core_jobs), run_unit))
  stopCluster(cl)
}
t_units <- (proc.time()[3] - t_start) / 60
cat(sprintf("units of work finished in %.1f minutes\n", t_units))
get <- function(id) readRDS(file.path(CACHE, paste0(id, ".rds")))

## ---------------------------------------------------------------------------
## 6. Summaries
## ---------------------------------------------------------------------------
within_of <- function(idx) if (length(idx)) mean(ss[idx]) else NA
jac <- function(a, b) { u <- length(union(a, b)); if (u) length(intersect(a, b)) / u else NA }
res_path <- function(name) file.path(OUT, "results", paste0("data_checks_", name, ".csv"))
core <- list(); for (v in c("orig", "garch")) for (w in 1:6)
  core[[v]][[wname[w]]] <- get(sprintf("core_%s_%s", v, wname[w]))
glas <- list(); mbr <- list(); for (v in c("orig", "garch")) for (w in 1:6) {
  glas[[v]][[wname[w]]] <- get(sprintf("glasso_%s_%s", v, wname[w]))
  mbr[[v]][[wname[w]]] <- get(sprintf("mb_%s_%s", v, wname[w]))
}

## 6.1 effective penalties (Table 6 with the corrected scale column), for the
##     original and the GARCH-standardized data, with GeneNet and beam
kap_tab <- do.call(rbind, lapply(c("orig", "garch"), function(v)
  do.call(rbind, lapply(1:6, function(w) {
    f <- core[[v]][[wname[w]]]
    data.frame(data = v, window = wname[w], n = f$n, lambda_cv = f$lambda_cv,
               kappa_hat = f$kappa_hat, s_hat = f$s_hat, c_hat = 1 / f$s_hat,
               kappa_floor_hat = f$kappa_floor_hat, sigma_max = f$sigma_max,
               kappa_over_sigma_max = f$kappa_hat / f$sigma_max,
               lambda_alt = f$lambda_alt,
               genenet_lambda = f$genenet_lambda, genenet_ridge = f$genenet_ridge,
               genenet_kappa = f$genenet_kappa, beam_alpha = f$beam_alpha,
               beam_ridge = f$beam_ridge, beam_kappa = f$beam_kappa,
               rt_kappa_corscale = f$rt_kappa_corscale,
               genenet_check = f$genenet_check, beam_check = f$beam_check)
  }))))
write.csv(kap_tab, res_path("kappa"), row.names = FALSE)

## 6.2 factor-adjusted effective penalties
fac_tab <- do.call(rbind, lapply(c("orig", "garch"), function(v)
  do.call(rbind, lapply(1:6, function(w) {
    f <- core[[v]][[wname[w]]]
    cbind(data = v, window = wname[w], lambda_cv = f$lambda_cv,
          kappa_hat = f$kappa_hat, f$factor)
  }))))
write.csv(fac_tab, res_path("factor_kappa"), row.names = FALSE)

## 6.3 leakage diagnostic per window. Plug-in rates use the extreme
##     eigenvalues of the full-sample covariance matrix of the same data,
##     with the smallest one corrected for its downward bias by the
##     parametric-bootstrap factor of data_extensions.R
lr_old <- read.csv(file.path(OUT, "results", "data_leakage_rates.csv"))
bias_min <- lr_old$bias_ratio_min[lr_old$adjustment == "none"][1]
bias_max <- lr_old$bias_ratio_max[lr_old$adjustment == "none"][1]
leak_tab <- do.call(rbind, lapply(c("orig", "garch"), function(v) {
  Xf <- make_Xv(v, 6); ev <- eigen(crossprod(Xf) / (nrow(Xf) - 1), TRUE,
                                   only.values = TRUE)$values
  smin_c <- min(ev) / bias_min; smax_c <- max(ev) / bias_max
  do.call(rbind, lapply(1:6, function(w) {
    f <- core[[v]][[wname[w]]]; k <- f$kappa_hat
    qr <- q_rates(k, smax_c, smin_c)
    rho <- (1 + k / smin_c) / (1 + k / smax_c)
    const <- smax_c * (1 + sqrt(rho))^2 / (2 * k * rho)
    ## smallest graph distance at which the bound (11) falls below 1/sqrt(n)
    dstar <- ceiling(log(1 / (sqrt(f$n) * const)) / log(qr["q_C"]))
    data.frame(data = v, window = wname[w], n = f$n, kappa_hat = k,
               kappa_floor_hat = f$kappa_floor_hat,
               kappa_over_sigma_max = k / f$sigma_max,
               sigma_min_corrected = smin_c, sigma_max_corrected = smax_c,
               shift_index = k / (k + smin_c),
               q_kappa = qr["q_kappa"], q_C = qr["q_C"], bound_const = const,
               bound_informative_from_distance = max(2, dstar),
               half_lambda = f$shift_half["lambda"], half_max = f$shift_half["max"],
               half_q99 = f$shift_half["q99"], half_q999 = f$shift_half["q999"],
               half_mean = f$shift_half["mean"], half_max_sqrt_n = sqrt(f$n) * f$shift_half["max"],
               half_jac100 = f$shift_half["jac100"], half_jac500 = f$shift_half["jac500"],
               half_spearman = f$shift_half["spearman"],
               double_max = f$shift_double["max"], double_q99 = f$shift_double["q99"],
               double_jac100 = f$shift_double["jac100"], double_jac500 = f$shift_double["jac500"],
               double_spearman = f$shift_double["spearman"], row.names = NULL)
  }))
}))
write.csv(leak_tab, res_path("leakage"), row.names = FALSE)

## 6.4 matched nominal against matched effective penalty: the full sample
##     refitted at the cross-validated penalty of each window, and at the
##     nominal penalty whose estimated effective penalty equals that of the
##     window. Window networks are those at their own cross-validated penalty.
mk_tab <- do.call(rbind, lapply(c("orig", "garch"), function(v) {
  Xf <- make_Xv(v, 6); nf <- nrow(Xf) - 1; spf <- spec_S(crossprod(Xf) / nf)
  fit_full <- function(l) {
    R <- pcor_of(ridge_Q(spf, l))
    list(lambda = l, kappa = kappa_hat(spf, l, nf), lfdr = lfdr_safe(R),
         top = rank_head(R, 1000), absr = abs(R[ut]))
  }
  own <- fit_full(core[[v]][["Full"]]$lambda_cv)
  do.call(rbind, lapply(1:5, function(w) {
    fw <- core[[v]][[wname[w]]]
    Xw <- make_Xv(v, w); nn <- nrow(Xw) - 1; spw <- spec_S(crossprod(Xw) / nn)
    Rw <- pcor_of(ridge_Q(spw, fw$lambda_cv)); rw <- abs(Rw[ut])
    win <- list(lfdr = fw$sel$RT_lFDR, top = fw$rank$RT[1:1000])
    fits <- list(own_cv = own, matched_lambda = fit_full(fw$lambda_cv),
                 matched_kappa = fit_full(lambda_for_kappa(spf, nf, fw$kappa_hat)))
    do.call(rbind, lapply(names(fits), function(nm) {
      ff <- fits[[nm]]
      data.frame(data = v, window = wname[w], comparison = nm,
                 window_lambda = fw$lambda_cv, window_kappa = fw$kappa_hat,
                 full_lambda = ff$lambda, full_kappa = ff$kappa,
                 window_lfdr_edges = length(win$lfdr), full_lfdr_edges = length(ff$lfdr),
                 window_lfdr_within = within_of(win$lfdr), full_lfdr_within = within_of(ff$lfdr),
                 jaccard_lfdr = jac(win$lfdr, ff$lfdr),
                 full_top100_within = within_of(ff$top[1:100]),
                 full_top500_within = within_of(ff$top[1:500]),
                 window_top500_within = within_of(win$top[1:500]),
                 jaccard_top100 = jac(win$top[1:100], ff$top[1:100]),
                 jaccard_top500 = jac(win$top[1:500], ff$top[1:500]),
                 spearman_abs_rho = cor(rw, ff$absr, method = "spearman"))
    }))
  }))
}))
write.csv(mk_tab, res_path("matched_kappa"), row.names = FALSE)

## 6.5 subsamples of several sizes: matched nominal against matched
##     effective penalty. The reference is the network of all rows at the
##     same nominal or the same effective penalty.
ref <- get(sprintf("sub_%d_%02d", nT, 1))
sub_raw <- do.call(rbind, lapply(NSUB[NSUB < nT], function(ns)
  do.call(rbind, lapply(1:20, function(r) {
    s <- get(sprintf("sub_%d_%02d", ns, r))
    do.call(rbind, lapply(names(s$at), function(a) {
      x <- s$at[[a]]; rf <- ref$at[[a]]
      if (is.null(x)) return(data.frame(nsub = ns, rep = r, match = a, floor = s$floor,
        lambda = NA, kappa = NA, lfdr_edges = NA, lfdr_within = NA, top500_within = NA,
        jaccard_lfdr_ref = NA, jaccard_top500_ref = NA, ref_kappa = rf$kappa))
      data.frame(nsub = ns, rep = r, match = a, floor = s$floor, lambda = x$lambda,
                 kappa = x$kappa, lfdr_edges = length(x$lfdr), lfdr_within = within_of(x$lfdr),
                 top500_within = within_of(x$top[1:500]),
                 jaccard_lfdr_ref = jac(x$lfdr, rf$lfdr),
                 jaccard_top500_ref = jac(x$top[1:500], rf$top[1:500]),
                 ref_kappa = rf$kappa)
    }))
  }))))
write.csv(sub_raw, res_path("subsample_reps"), row.names = FALSE)
agg <- function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x))))
sub_tab <- do.call(rbind, lapply(split(sub_raw, list(sub_raw$match, sub_raw$nsub), drop = TRUE),
  function(d) {
    vars <- c("floor", "lambda", "kappa", "lfdr_edges", "lfdr_within", "top500_within",
              "jaccard_lfdr_ref", "jaccard_top500_ref")
    out <- data.frame(nsub = d$nsub[1], match = d$match[1], reps_ok = sum(!is.na(d$kappa)))
    for (vv in vars) { a <- agg(d[[vv]]); out[[vv]] <- a["mean"]; out[[paste0(vv, "_se")]] <- a["se"] }
    out
  }))
ref_rows <- do.call(rbind, lapply(names(ref$at), function(a) {
  x <- ref$at[[a]]
  data.frame(nsub = nT, match = a, reps_ok = 1, floor = ref$floor, lambda = x$lambda,
             kappa = x$kappa, lfdr_edges = length(x$lfdr), lfdr_within = within_of(x$lfdr),
             top500_within = within_of(x$top[1:500]))
}))
sub_tab <- merge(sub_tab, ref_rows, all = TRUE)
sub_tab <- sub_tab[order(sub_tab$match, sub_tab$nsub), ]
write.csv(sub_tab, res_path("subsample"), row.names = FALSE)

## 6.6 tuning of the graphical lasso and of neighbourhood selection
tune_rows <- list()
add_tune <- function(v, method, tuning, getter, info = function(w) c(NA, NA)) {
  el <- lapply(1:6, getter)
  ne <- sapply(el[1:5], length)
  inf <- sapply(1:6, function(w) { x <- info(w); if (length(x) < 2) c(NA, x) else x })
  tune_rows[[length(tune_rows) + 1]] <<- data.frame(data = v, method = method,
    tuning = tuning, W1 = ne[1], W2 = ne[2], W3 = ne[3], W4 = ne[4], W5 = ne[5],
    edges_mean = mean(ne), within_sector = mean(sapply(el[1:5], within_of), na.rm = TRUE),
    jaccard_consecutive = mean(sapply(1:4, function(w) jac(el[[w]], el[[w + 1]])), na.rm = TRUE),
    edges_full = length(el[[6]]), within_full = within_of(el[[6]]),
    opt_index_windows = paste(inf[1, 1:5], collapse = " "),
    opt_index_full = inf[1, 6], path_length = inf[2, 6])
}
sel_main <- readRDS(file.path(OUT, "results", "data_selections.rds"))
for (v in c("orig", "garch")) {
  G <- glas[[v]]; M <- mbr[[v]]
  if (v == "orig")
    add_tune(v, "GLASSO", "EBIC, main path (30, 0.05)", function(w)
      if (w < 6) edges_of(sel_main[[wname[w]]]$GLASSO) else G[[w]]$old$sel,
      function(w) c(if (w < 6) NA else G[[w]]$old$opt, 30))
  add_tune(v, "GLASSO", "EBIC, finer path (60, 0.1 windows / 0.02 full)",
           function(w) G[[w]]$fine$sel, function(w) c(G[[w]]$fine$opt, length(G[[w]]$fine$lambda)))
  add_tune(v, "GLASSO", "StARS (beta 0.1, path 30, 0.2)", function(w) G[[w]]$stars$sel,
           function(w) c(G[[w]]$stars$opt, length(G[[w]]$stars$lambda)))
  add_tune(v, "MB", "RIC, main path (30, 0.05)", function(w) M[[w]]$old$sel,
           function(w) c(M[[w]]$old$opt, length(M[[w]]$old$lambda)))
  add_tune(v, "MB", "RIC, finer path (60, 0.01)", function(w) M[[w]]$ric_fine$sel,
           function(w) c(M[[w]]$ric_fine$opt, 60))
  add_tune(v, "MB", "StARS (beta 0.1, path 60, 0.01)", function(w) M[[w]]$stars$sel,
           function(w) c(M[[w]]$stars$opt, 60))
  add_tune(v, "MB", "nodewise EBIC (gamma 0.5), OR rule", function(w) M[[w]]$ebic$or)
  add_tune(v, "MB", "nodewise EBIC (gamma 0.5), AND rule", function(w) M[[w]]$ebic$and)
}
tune_tab <- do.call(rbind, tune_rows)
write.csv(tune_tab, res_path("tuning"), row.names = FALSE)

## 6.7 the network comparison (Table 7) with GeneNet, beam and the retuned
##     lasso methods, for the original and the GARCH-standardized data
meth_rows <- list()
for (v in c("orig", "garch")) {
  for (m in names(core[[v]][[1]]$sel)) {
    el <- lapply(1:6, function(w) core[[v]][[wname[w]]]$sel[[m]])
    ne <- sapply(el[1:5], length)
    meth_rows[[length(meth_rows) + 1]] <- data.frame(data = v, method = m,
      edges_mean = mean(ne), within_sector = mean(sapply(el[1:5], within_of), na.rm = TRUE),
      jaccard_consecutive = mean(sapply(1:4, function(w) jac(el[[w]], el[[w + 1]])), na.rm = TRUE),
      edges_full = length(el[[6]]), within_full = within_of(el[[6]]))
  }
}
meth_tab <- rbind(do.call(rbind, meth_rows),
                  tune_tab[, c("data", "method", "edges_mean", "within_sector",
                               "jaccard_consecutive", "edges_full", "within_full")])
meth_tab$method[(length(meth_rows) + 1):nrow(meth_tab)] <-
  paste0(tune_tab$method, ": ", tune_tab$tuning)
write.csv(meth_tab, res_path("methods"), row.names = FALSE)

## 6.8 reproducibility at matched sizes: every ranking cut at K edges in each
##     window, with K fixed or K equal to the RT-MPS count of the window; the
##     lasso paths are cut at the member closest to K (finer paths)
closest <- function(sizes, K) which.min(abs(sizes - K))
ms_rows <- list()
for (v in c("orig", "garch")) {
  Kmps <- sapply(1:5, function(w) length(core[[v]][[wname[w]]]$sel$RT_MPS))
  for (Kname in c("RT-MPS size", KS)) {
    Kw <- if (Kname == "RT-MPS size") Kmps else rep(as.numeric(Kname), 5)
    sets <- list()
    for (m in names(core[[v]][[1]]$rank))
      sets[[m]] <- lapply(1:5, function(w) core[[v]][[wname[w]]]$rank[[m]][seq_len(Kw[w])])
    sets$GLASSO <- lapply(1:5, function(w) { g <- glas[[v]][[wname[w]]]$fine
      g$path[[closest(g$size, Kw[w])]] })
    sets$MB <- lapply(1:5, function(w) { g <- mbr[[v]][[wname[w]]]$fine
      g$path[[closest(g$size, Kw[w])]] })
    for (m in names(sets)) ms_rows[[length(ms_rows) + 1]] <- data.frame(data = v,
      K = Kname, method = m, mean_size = mean(sapply(sets[[m]], length)),
      within_sector = mean(sapply(sets[[m]], within_of)),
      jaccard_consecutive = mean(sapply(1:4, function(w) jac(sets[[m]][[w]], sets[[m]][[w + 1]]))))
  }
}
write.csv(do.call(rbind, ms_rows), res_path("matched_size"), row.names = FALSE)

## 6.9 semi-synthetic check: does rho at the estimated effective penalty
##     centre the ridge partial correlations under a spiked Sigma, with
##     Gaussian and with t rows?
Xf <- make_Xv("orig", 6); Sig <- cov2cor(crossprod(Xf))
eS <- eigen(Sig, symmetric = TRUE)
rho_pop <- function(k) pcor_of(eS$vectors %*% (t(eS$vectors) / (eS$values + k)))[ut]
P_theta <- pcor_of(solve(Sig))[ut]
sim_rows <- list()
for (d in c("gauss", "t5", "t5ns")) {
  ch <- lapply(1:4, function(c) get(sprintf("sim_%s_%d", d, c)))
  B <- sum(sapply(ch, `[[`, "reps"))
  kap <- do.call(rbind, lapply(ch, `[[`, "kappa"))
  flo <- unlist(lapply(ch, `[[`, "floor"))
  shift <- do.call(rbind, lapply(ch, `[[`, "shift"))
  for (g in seq_along(SIM_LAM)) {
    mu <- Reduce(`+`, lapply(ch, function(c) c$acc[[g]]$sum)) / B
    sq <- Reduce(`+`, lapply(ch, function(c) c$acc[[g]]$sq)) / B
    sdv <- sqrt(pmax(sq - mu^2, 0))
    rms <- function(k) sqrt(mean((mu - rho_pop(k))^2))
    opt <- optimize(function(lk) rms(exp(lk)), log(c(1e-3, 1e3)))
    k_fit <- exp(opt$minimum); kh <- mean(kap[, g])
    pk <- rho_pop(kh); ph <- rho_pop(kh / 2)
    sim_rows[[length(sim_rows) + 1]] <- data.frame(rows = d, reps = B, lambda = SIM_LAM[g],
      kappa_hat_mean = kh, kappa_hat_sd = sd(kap[, g]), kappa_fit = k_fit,
      ratio_hat_fit = kh / k_fit, rms_fit = opt$objective, rms_at_kappa_hat = rms(kh),
      rms_sd_rho = sqrt(mean(sdv^2)), rms_at_zero = sqrt(mean((mu - P_theta)^2)),
      floor_mean = mean(flo),
      pop_half_max = if (g == 2) max(abs(pk - ph)) else NA,
      sample_half_max = if (g == 2) mean(shift[, 1]) else NA,
      pop_total_max = if (g == 2) max(abs(pk - P_theta)) else NA,
      pop_total_q99 = if (g == 2) unname(quantile(abs(pk - P_theta), 0.99)) else NA,
      pop_half_q99 = if (g == 2) unname(quantile(abs(pk - ph), 0.99)) else NA,
      sample_half_q99 = if (g == 2) mean(shift[, 2]) else NA)
  }
}
write.csv(do.call(rbind, sim_rows), res_path("semisynthetic"), row.names = FALSE)

## 6.10 GARCH fits and running times
write.csv(data.frame(stock = stockdata$info[, 1], sector = sector, method = gres$method,
                     gres$coef, row.names = NULL), res_path("garch_fits"), row.names = FALSE)
rt <- do.call(rbind, lapply(c(names(jobs), sapply(core_jobs, `[[`, "id")), function(id)
  data.frame(unit = id, minutes = get(id)$seconds / 60)))
write.csv(rt, res_path("runtime"), row.names = FALSE)
writeLines(c(sprintf("Run time of this call: %.1f minutes on %d cores (units of work: %.1f minutes)",
                     (proc.time()[3] - t_start) / 60, NCORES, t_units),
             sprintf("Sum of the run times of all units of work: %.1f minutes", sum(rt$minutes)),
             capture.output(sessionInfo())),
           file.path(OUT, "results", "sessionInfo_checks.txt"))
cat("Further data analyses finished.\n")
