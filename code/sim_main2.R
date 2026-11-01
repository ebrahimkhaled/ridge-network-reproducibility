###############################################################################
## Second extension of the main simulation design for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Purpose: main design (p = 100; band, random, hub, scale-free graphs;
##   n = 50, 100, 200, 400), extended to 500 replicates.
##   Replicates 1-100 use the random numbers of sim_ridge_support.R: the data
##   and the cross-validation folds are rebuilt from the stream states stored
##   by sim_extensions.R in ../results/ext_work/, and the published
##   quantities are read from there. Replicates 101-500 run the full
##   one_rep_ext() of sim_extensions.R on new streams.
##   Added in every replicate:
##     - competitors: GeneNet (shrinkage partial correlations, local fdr at
##       0.8 and 0.9, Benjamini-Hochberg and tail-area Fdr at 0.1), beam
##       (Leday and Richardson, 2019, default settings), the vWP estimator
##       with the default target of rags2ridges, and (replicates 1-100 only)
##       the ridge and PLS networks of parcor (Kraemer et al., 2009);
##     - every ranking at the true number of edges and at the size selected
##       by Liu (2013); FDR rules at the nominal level 0.1;
##     - the plain ridge ranking at a 100-point CV grid, at multiples of the
##       CV penalty and at a subsampling-stability penalty;
##     - the oracle ranking with the population leakage removed;
##     - the extrapolated ranking with lambda_2 from a fixed grid and at
##       large first penalties, with its population target;
##     - signed partial correlations (replicates 1-100) for the fit of the
##       population ridge network to other ridge-type estimators.
##   A second part rebuilds the p = 400 block of sim_extensions.R (same
##   random numbers) and adds the extrapolation, the oracle ranking and the
##   cross-validation curve.
##
## Requirements: R >= 4.1; packages huge, glmnet, fdrtool, igraph, GGMridge,
##   SILGGM, GeneNet, beam, rags2ridges, parcor (archived on CRAN, version
##   0.2-6). An extra library folder can be given in RIDGE_RLIB.
## Usage  : Rscript sim_main2.R [main] [hd]      (default: both), from code/.
##   Every replicate is written to ../results/ext_work3/ when finished; a
##   restarted run skips the replicates already there. RIDGE_CORES sets the
##   number of worker processes (default 11). sim_main2_summary.R makes the
##   tables.
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("main", "hd")
OUT   <- ".."
WORK0 <- normalizePath(file.path(OUT, "results", "ext_work"))
dir.create(file.path(OUT, "results", "ext_work3"), showWarnings = FALSE)
WORK  <- normalizePath(file.path(OUT, "results", "ext_work3"))
NCORE <- min(as.integer(Sys.getenv("RIDGE_CORES", "11")), detectCores() - 1)
CODE  <- normalizePath(".")
NREP  <- 500       # replicates per setting
NOLD  <- 100       # replicates with the published random numbers
NPAR  <- 100       # replicates in which parcor is run (it is slow)

## Packages and the function definitions of the two earlier scripts (only the
## definitions are evaluated, the studies are not run).
load_tools2 <- function(code_dir, lib_extra) {
  Sys.setenv(OMP_NUM_THREADS = 1)
  if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
  suppressPackageStartupMessages({
    library(huge); library(glmnet); library(fdrtool); library(igraph)
    library(GGMridge); library(SILGGM); library(GeneNet); library(beam)
    library(parcor)
  })
  for (f in c("sim_ridge_support.R", "sim_extensions.R")) {
    for (e in parse(file.path(code_dir, f))) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function")))
        eval(e, envir = .GlobalEnv)
    }
  }
  assign("MULT", c(0.3, 1, 3, 10), envir = .GlobalEnv)
  assign("NPERM", 99, envir = .GlobalEnv)
  assign("ALPHA", 0.1, envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG")
  invisible(TRUE)
}
load_tools2(CODE, lib_extra)

p_sim  <- 100
n_grid <- c(50, 100, 200, 400)
graphs <- c("band", "random", "hub", "scale-free")
GRID30  <- exp(seq(log(1e-2), log(20), length.out = 30))
GRID100 <- exp(seq(log(1e-2), log(20), length.out = 100))
## wider grid for the vWP estimator with the default target: the target can
## be large (mean inverse eigenvalue) and the CV penalty then small
GRIDT   <- exp(seq(log(1e-5), log(20), length.out = 50))

## ---------------------------------------------------------------------------
## 1. Helpers
## ---------------------------------------------------------------------------

## stream i of a seed: nextRNGStream() applied i times
stream_of <- function(seed, i) {
  set.seed(seed)
  s <- .Random.seed
  for (k in seq_len(i)) s <- parallel::nextRNGStream(s)
  s
}

## Cross-validated loss curve with given folds. Same arithmetic as
## cv_lambda() of sim_ridge_support.R; "alt" is the vWP estimator with target
## tgt * I (tgt = 0: null target):
## [ {lambda I + (S - lambda T)^2 / 4}^{1/2} + (S - lambda T) / 2 ]^{-1}.
cv_curve <- function(X, fold, grid, type = "ridge", tgt = 0) {
  K <- max(fold)
  loss <- matrix(0, K, length(grid))
  for (k in 1:K) {
    Xtr <- X[fold != k, , drop = FALSE]
    Xte <- X[fold == k, , drop = FALSE]
    sp <- spec_S(crossprod(Xtr) / nrow(Xtr))
    Ste <- crossprod(Xte) / nrow(Xte)
    ste_diag <- colSums(sp$vec * (Ste %*% sp$vec))
    for (g in seq_along(grid)) {
      dv <- if (type == "ridge") 1 / (sp$val + grid[g])
            else { a <- sp$val - grid[g] * tgt
                   1 / (sqrt(grid[g] + a^2 / 4) + a / 2) }
      loss[k, g] <- sum(ste_diag * dv) - sum(log(dv))
    }
  }
  colMeans(loss)
}

## vWP estimator with target tgt * I, from the spectrum of S
altridge_T <- function(sp, lambda, tgt) {
  a <- sp$val - lambda * tgt
  d <- sqrt(lambda + a^2 / 4) + a / 2
  sp$vec %*% (t(sp$vec) / d)
}

## default target of rags2ridges, default.target(S) with type "DAIE":
## mean of the inverse eigenvalues not below 1e-4 times the largest
daie_target <- function(sp) {
  ev <- sp$val
  mean(1 / ev[ev >= max(ev) * 1e-4])
}

## local fdr fit of fdrtool on the upper triangle of a partial-correlation
## matrix (the same call as lfdr_select())
lfdr_fit <- function(R) {
  ut <- upper.tri(R)
  f <- fdrtool(R[ut], statistic = "correlation", plot = FALSE, verbose = FALSE)
  list(lfdr = f$lfdr, qval = f$qval, pval = f$pval,
       eta0 = unname(f$param[1, "eta0"]))
}
vec_to_adj <- function(v, p) {
  A <- matrix(FALSE, p, p); A[upper.tri(A)] <- v
  A | t(A)
}

## Ranking summaries of a score vector s (upper triangle, larger = stronger):
## best MCC over all thresholds, MCC and TPR at the m_true strongest pairs,
## MCC, TPR and FDP at the m_liu strongest pairs, AUC.
rank_summary <- function(s, e, m_true, m_liu) {
  o <- order(s, decreasing = TRUE)
  tp <- cumsum(e[o]); fp <- cumsum(!e[o])
  P <- sum(e); N <- sum(!e); fn <- P - tp; tn <- N - fp
  den <- sqrt(tp + fp) * sqrt(tp + fn) * sqrt(tn + fp) * sqrt(tn + fn)
  mcc <- ifelse(den > 0, (tp * tn - fp * fn) / den, 0)
  at <- function(k) if (k < 1) c(NA, NA, NA) else
    c(mcc[k], tp[k] / P, fp[k] / k)
  auc <- (sum(rank(s)[e]) - P * (P + 1) / 2) / (P * N)
  a1 <- at(m_true); a2 <- at(m_liu)
  c(best = max(mcc), mcc_m = a1[1], tpr_m = a1[2], mcc_liu = a2[1],
    tpr_liu = a2[2], fdp_liu = a2[3], auc = auc)
}

## The same for a path of graphs: the graph whose size is closest to the
## target size (the smaller one when two are equally close).
path_summary <- function(path, E, m_true, m_liu) {
  ut <- upper.tri(E)
  sz <- sapply(path, function(A) sum(A[ut]))
  met <- t(sapply(path, function(A) { diag(A) <- FALSE; sel_metrics(A, E) }))
  pick <- function(k) {
    if (k < 1) return(c(NA, NA, NA, NA))
    d <- abs(sz - k); j <- which(d == min(d)); j <- j[which.min(sz[j])]
    c(met[j, "MCC"], met[j, "TPR"], met[j, "FDP"], sz[j])
  }
  a1 <- pick(m_true); a2 <- pick(m_liu)
  c(best = max(met[, "MCC"]), mcc_m = a1[1], tpr_m = a1[2], size_m = a1[4],
    mcc_liu = a2[1], tpr_liu = a2[2], fdp_liu = a2[3], size_liu = a2[4])
}

## margin and best MCC of a population partial-correlation matrix
pop_summary <- function(R, E) {
  m <- pop_margin(R, E)
  c(zeta = unname(m["zeta"]), ell = unname(m["ell"]),
    mcc = oracle_mcc_threshold(R, E))
}

## Two-point extrapolation with a given second penalty (Section 4.2 of the
## paper); returns the estimate and the ratio h = kappa_2 / kappa_1.
richardson_Q2 <- function(sp, n, lambda1, lambda2) {
  k1 <- kappa_hat(sp, lambda1, n); k2 <- kappa_hat(sp, lambda2, n)
  w1 <- k2 / (k2 - k1); w2 <- -k1 / (k2 - k1)
  list(T = w1 * s_hat(sp, lambda1, n) * ridge_Q(sp, lambda1) +
           w2 * s_hat(sp, lambda2, n) * ridge_Q(sp, lambda2),
       kappa = c(k1, k2), h = k2 / k1)
}
## second penalty from a fixed grid: the grid point above lambda1 whose
## kappa_hat is closest to 2 kappa_hat(lambda1)
grid_lambda2 <- function(sp, n, lambda1, grid) {
  g <- grid[grid > lambda1 * (1 + 1e-8)]
  if (!length(g)) return(NA)
  k1 <- kappa_hat(sp, lambda1, n)
  kk <- sapply(g, kappa_hat, sp = sp, n = n)
  g[which.min(abs(kk - 2 * k1))]
}
## population target of the extrapolation at effective penalties k1 < k2
pop_extrap <- function(Sig, k1, k2) {
  p <- ncol(Sig)
  w1 <- k2 / (k2 - k1); w2 <- -k1 / (k2 - k1)
  w1 * solve(Sig + k1 * diag(p)) + w2 * solve(Sig + k2 * diag(p))
}

## Subsampling stability of the top-K pairs of the ridge ranking. For each
## penalty omega * tr(S) / p, B subsamples of size n / 2 are drawn and the
## mean pairwise Jaccard index of their K strongest pairs is computed; the
## penalty with the largest stability is chosen.
stab_select <- function(X, omegas, Ks, B = 20) {
  n <- nrow(X); p <- ncol(X); ut <- upper.tri(diag(p))
  lams <- omegas * sum(X^2) / (n - 1) / p
  tops <- lapply(Ks, function(K)
    replicate(length(lams), vector("list", B), simplify = FALSE))
  for (b in 1:B) {
    idx <- sample(n, floor(n / 2))
    Xb <- scale(X[idx, , drop = FALSE], center = TRUE, scale = FALSE)
    spb <- spec_S(crossprod(Xb) / (nrow(Xb) - 1))
    for (g in seq_along(lams)) {
      s <- abs(pcor_of(ridge_Q(spb, lams[g]))[ut])
      o <- order(s, decreasing = TRUE)
      for (k in seq_along(Ks)) tops[[k]][[g]][[b]] <- o[seq_len(Ks[k])]
    }
  }
  jac <- function(L) {
    v <- c()
    for (a in 1:(B - 1)) for (b in (a + 1):B)
      v <- c(v, length(intersect(L[[a]], L[[b]])) / length(union(L[[a]], L[[b]])))
    mean(v)
  }
  st <- sapply(seq_along(Ks), function(k) sapply(tops[[k]], jac))
  st <- matrix(st, nrow = length(lams))
  list(stab = st, lambda = lams,
       pick = apply(st, 2, function(v) lams[which.max(v)]))
}

## ---------------------------------------------------------------------------
## 2. Main design: one replicate
## ---------------------------------------------------------------------------
main2_rep <- function(job) {
  g <- job$graph; n <- job$n; r <- job$rep
  f_out <- file.path(job$work, sprintf("main2_%s_%d_rep%03d.rds", g, n, r))
  if (file.exists(f_out)) return("skip")
  t0 <- proc.time()[3]
  gi <- match(g, graphs); ni <- match(n, n_grid)
  mod <- make_model(g, p_sim)            # uses set.seed() internally
  Sig <- mod$Sig; Theta <- mod$Theta; E <- mod$E; p <- ncol(Sig)
  ut <- upper.tri(E); e <- E[ut]; m_true <- sum(e)

  ## stream state before the replicate
  if (r <= NOLD) {
    s0 <- if (r <= 2) stream_of(1000 + 10 * gi + ni, r)
          else readRDS(file.path(job$work0, sprintf("main_%s_%d_rep%03d.rds",
                                                    g, n, r - 2)))$seed_after
    ext <- readRDS(file.path(job$work0, sprintf("main_%s_%d_rep%03d.rds",
                                                g, n, r)))$res
    ext_new <- NULL
  } else {
    s0 <- stream_of(5000 + 10 * gi + ni, r - NOLD)
    assign(".Random.seed", s0, envir = .GlobalEnv)
    ext <- one_rep_ext(n, Sig, Theta, E, mod$D, job$lam_fix)
    ext$rho_hat <- ext$rho_kap <- NULL
    ext_new <- ext
  }

  ## rebuild the data, the CV folds and the CV penalties of the replicate
  assign(".Random.seed", s0, envir = .GlobalEnv)
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sig)
  X <- scale(X, center = TRUE, scale = FALSE)
  neff <- n - 1
  S <- crossprod(X) / neff
  sp <- spec_S(S)
  s_f1 <- .Random.seed
  lam <- cv_lambda(X, GRID30, type = "ridge")
  s_f2 <- .Random.seed
  lam_a <- cv_lambda(X, GRID30, type = "alt")
  assign(".Random.seed", s_f1, envir = .GlobalEnv)
  fold1 <- sample(rep_len(1:5, n))
  assign(".Random.seed", s_f2, envir = .GlobalEnv)
  fold2 <- sample(rep_len(1:5, n))
  ## all further random numbers come from an auxiliary stream
  assign(".Random.seed", stream_of(7000 + 10 * gi + ni, r), envir = .GlobalEnv)

  R_plain <- pcor_of(ridge_Q(sp, lam))
  kap <- kappa_hat(sp, lam, neff)
  RR <- richardson_Q(sp, neff, lam, cfac = 2)
  R_rich <- pcor_of(RR$T)
  R_alt <- pcor_of(altridge_Q(sp, lam_a))
  R_kap <- pcor_of(solve(Sig + kap * diag(p)))
  P_true <- pcor_of(Theta)
  C <- cov2cor(S)
  spC <- spec_S(C)

  ## checks against the stored replicate
  cv1 <- cv_curve(X, fold1, GRID30, "ridge")
  chk <- c(lambda = lam - ext$info[["lambda"]],
           kappa = kap - ext$info[["kappa"]],
           lam_a = abs(GRID30[which.min(cv_curve(X, fold2, GRID30, "alt"))] - lam_a),
           fold1 = abs(GRID30[which.min(cv1)] - lam),
           orc_RT = oracle_mcc_threshold(R_plain, E) - ext$orc[["RT"]],
           orc_vWP = oracle_mcc_threshold(R_alt, E) - ext$orc[["vWP"]])

  ## ---- tuning variants of the plain ridge ranking ----
  cv100 <- cv_curve(X, fold1, GRID100, "ridge")
  lam100 <- GRID100[which.min(cv100)]
  R_cv100 <- pcor_of(ridge_Q(sp, lam100))
  lf_rt <- lfdr_fit(R_plain)
  K_eta <- min(1000, max(20, round((1 - lf_rt$eta0) * sum(ut))))
  omegas <- c(0.1, 0.3, 1, 3, 10, 30)
  st <- stab_select(X, omegas, Ks = c(K_eta, p))
  R_stab <- pcor_of(ridge_Q(sp, st$pick[1]))
  R_stabp <- pcor_of(ridge_Q(sp, st$pick[2]))
  R_x <- lapply(c(3, 10, 30), function(cc) pcor_of(ridge_Q(sp, cc * lam)))

  ## ---- vWP with the default target of rags2ridges ----
  tgt <- daie_target(sp)
  cvt <- cv_curve(X, fold2, GRIDT, "alt", tgt = tgt)
  lam_t <- GRIDT[which.min(cvt)]
  R_altT <- pcor_of(altridge_T(sp, lam_t, tgt))

  ## ---- GeneNet: Schafer-Strimmer shrinkage partial correlations ----
  pc_gn <- ggm.estimate.pcor(X, method = "static", verbose = FALSE)
  lam_gn <- attr(pc_gn, "lambda")
  R_gn <- matrix(as.numeric(pc_gn), p, p)
  invisible(capture.output(te <- network.test.edges(pc_gn, plot = FALSE,
                                                     verbose = FALSE)))
  idx <- cbind(pmin(te$node1, te$node2), pmax(te$node1, te$node2))
  gn_prob <- gn_q <- gn_p <- matrix(NA, p, p)
  gn_prob[idx] <- te$prob; gn_q[idx] <- te$qval; gn_p[idx] <- te$pval
  gn_prob <- gn_prob[ut]; gn_q <- gn_q[ut]; gn_p <- gn_p[ut]
  ## (1 - l) R + l I = (1 - l) {R + l / (1 - l) I}: ridge on R at l / (1 - l)
  kap_gn <- if (lam_gn < 1) kappa_hat(spC, lam_gn / (1 - lam_gn), neff) else Inf

  ## ---- beam (default settings) ----
  bf <- beam(X, type = "conditional", verbose = FALSE)
  btab <- bf@table
  bsel <- beam.select(bf, thres = 0.1, method = "BH")
  bm <- as.matrix(bsel@conditional)
  b_bh <- p.adjust(btab[, "p_tail_prob"], method = "BH") < 0.1
  ## beam's partial correlations are ridge partial correlations of the
  ## correlation matrix at penalty (delta - p - 1) / n
  lam_b <- (bf@deltaOpt - p - 1) / n
  b_chk <- max(abs(pcor_of(solve(cor(X) + lam_b * diag(p)))[ut] - btab[, "p_cor"]))
  kap_b <- kappa_hat(spC, lam_b, neff)
  ## ranking by the tail probability, ties broken by |partial correlation|
  b_score <- numeric(sum(ut))
  b_score[order(btab[, "p_tail_prob"], -abs(btab[, "p_cor"]))] <- rev(seq_len(sum(ut)))

  ## ---- Ha-Sun ranking (ridge of the correlation matrix at its penalty) ----
  lam_hs <- ext$hs[["lambda"]]
  R_hs <- pcor_of(solve(C + lam_hs * diag(p)))

  ## ---- parcor (replicates 1-100) ----
  par_ok <- r <= NPAR
  if (par_ok) {
    invisible(capture.output(pr <- tryCatch(ridge.net(X, k = 10)$pcor,
                                            error = function(e) NULL)))
    invisible(capture.output(pl <- tryCatch(pls.net(X, k = 10, ncomp = 15)$pcor,
                                            error = function(e) NULL)))
  }

  ## ---- oracle rankings with the population leakage removed ----
  R_orc <- R_plain - (R_kap - P_true)       # leakage and shrinkage removed
  R_orl <- R_plain; R_orl[!E] <- R_plain[!E] - R_kap[!E]   # leakage only
  diag(R_orl) <- 1

  ## ---- extrapolation: second penalty from fixed grids ----
  l2_30 <- grid_lambda2(sp, neff, lam, GRID30)
  l2_100 <- grid_lambda2(sp, neff, lam, GRID100)
  RG30 <- richardson_Q2(sp, neff, lam, l2_30)
  RG100 <- richardson_Q2(sp, neff, lam, l2_100)
  l2_root <- lambda_for_kappa(sp, neff, 2 * kap, lower = lam)
  re_info <- c(lambda1 = lam, lambda2 = l2_root, ratio = l2_root / lam,
               pd_guarantee = l2_root <= 4 * lam,
               min_eig = min(eigen(RR$T, TRUE, only.values = TRUE)$values),
               lambda2_g30 = l2_30, h_g30 = RG30$h, lambda2_g100 = l2_100,
               h_g100 = RG100$h,
               pd_g30 = min(eigen(RG30$T, TRUE, only.values = TRUE)$values) > 0,
               kappa1 = RR$kappa[1], kappa2 = RR$kappa[2])
  ## population targets of RT and RE at lambda_CV and at large penalties
  pop_re <- list()
  for (om in c(0, 1, 3, 10)) {
    l1 <- if (om == 0) lam else om * sum(diag(S)) / p
    k1 <- kappa_hat(sp, l1, neff)
    l2 <- lambda_for_kappa(sp, neff, 2 * k1, lower = l1)
    k2 <- kappa_hat(sp, l2, neff)
    Rre <- pcor_of(richardson_Q(sp, neff, l1, cfac = 2)$T)
    pop_re[[paste0("om", om)]] <- c(
      lambda1 = l1, kappa1 = k1, ratio = l2 / l1,
      RT_best = oracle_mcc_threshold(pcor_of(ridge_Q(sp, l1)), E),
      RE_best = oracle_mcc_threshold(Rre, E),
      setNames(pop_summary(pcor_of(solve(Sig + k1 * diag(p))), E),
               c("zeta_RT", "ell_RT", "mcc_RT")),
      setNames(pop_summary(pcor_of(pop_extrap(Sig, k1, k2)), E),
               c("zeta_RE", "ell_RE", "mcc_RE")))
  }
  pop_re <- do.call(rbind, pop_re)

  ## ---- rankings ----
  m_liu <- ext$met2["Liu2013", "NE"]
  scores <- list(RT = abs(R_plain[ut]), RE = abs(R_rich[ut]),
                 vWP = abs(R_alt[ut]), vWPt = abs(R_altT[ut]),
                 GeneNet = abs(R_gn[ut]), beam = b_score, HS = abs(R_hs[ut]),
                 Marginal = abs(C[ut]), RT_cv100 = abs(R_cv100[ut]),
                 RT_stab = abs(R_stab[ut]), RT_stabp = abs(R_stabp[ut]),
                 RT_x3 = abs(R_x[[1]][ut]), RT_x10 = abs(R_x[[2]][ut]),
                 RT_x30 = abs(R_x[[3]][ut]),
                 RE_g30 = abs(pcor_of(RG30$T)[ut]),
                 RE_g100 = abs(pcor_of(RG100$T)[ut]),
                 ORC = abs(R_orc[ut]), ORC_leak = abs(R_orl[ut]))
  if (par_ok && !is.null(pr)) scores$parcor_ridge <- abs(pr[ut])
  if (par_ok && !is.null(pl)) scores$parcor_pls <- abs(pl[ut])
  rank <- t(sapply(scores, rank_summary, e = e, m_true = m_true, m_liu = m_liu))

  gl <- glasso_ebic(X, S)
  hmb <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.01,
              verbose = FALSE)
  mbp <- lapply(hmb$path, function(A) as.matrix(A) != 0)
  paths <- rbind(GLASSO = path_summary(gl$path, E, m_true, m_liu),
                 MB = path_summary(mbp, E, m_true, m_liu))
  chk <- c(chk, orc_GL = paths["GLASSO", "best"] - ext$orc[["GLASSO"]],
           orc_MB = paths["MB", "best"] - ext$orc[["MB"]])

  ## ---- new selections ----
  lf_a <- lfdr_fit(R_alt); lf_t <- lfdr_fit(R_altT)
  sel <- list(
    vWP_lFDR_chk = vec_to_adj(1 - lf_a$lfdr >= 0.9, p),
    vWPt_lFDR90  = vec_to_adj(1 - lf_t$lfdr >= 0.9, p),
    vWPt_lFDR99  = vec_to_adj(1 - lf_t$lfdr >= 0.99, p),
    RT_Fdr10     = vec_to_adj(lf_rt$qval < 0.1, p),
    vWP_Fdr10    = vec_to_adj(lf_a$qval < 0.1, p),
    vWPt_Fdr10   = vec_to_adj(lf_t$qval < 0.1, p),
    RT_BH10      = vec_to_adj(p.adjust(lf_rt$pval, "BH") < 0.1, p),
    GeneNet_lFDR80 = vec_to_adj(gn_prob >= 0.8, p),
    GeneNet_lFDR90 = vec_to_adj(gn_prob >= 0.9, p),
    GeneNet_BH10 = vec_to_adj(p.adjust(gn_p, "BH") < 0.1, p),
    GeneNet_Fdr10 = vec_to_adj(gn_q < 0.1, p),
    beam_BH10    = vec_to_adj(b_bh, p),
    RT_cv100_lFDR90 = lfdr_select(R_cv100),
    RT_stab_lFDR90  = lfdr_select(R_stab))
  if (par_ok && !is.null(pr)) {
    lf <- lfdr_fit(pr)
    sel$parcor_ridge_lFDR80 <- vec_to_adj(1 - lf$lfdr >= 0.8, p)
    sel$parcor_ridge_lFDR90 <- vec_to_adj(1 - lf$lfdr >= 0.9, p)
  }
  if (par_ok && !is.null(pl)) {
    lf <- lfdr_fit(pl)
    sel$parcor_pls_lFDR80 <- vec_to_adj(1 - lf$lfdr >= 0.8, p)
    sel$parcor_pls_lFDR90 <- vec_to_adj(1 - lf$lfdr >= 0.9, p)
  }
  for (k in names(sel)) diag(sel[[k]]) <- FALSE
  met3 <- t(sapply(sel, sel_metrics, E = E))
  chk <- c(chk, sel_vWP = sum(sel$vWP_lFDR_chk[ut]) - ext$met["vWP_lFDR", "NE"],
           beam_ridge = b_chk, beam_nsel = nrow(bm) - sum(b_bh))

  tune <- c(lambda = lam, kappa = kap, lambda_cv100 = lam100,
            kappa_cv100 = kappa_hat(sp, lam100, neff),
            lambda_vWP = lam_a, lambda_vWPt = lam_t, target_vWPt = tgt,
            vWPt_at_edge = lam_t %in% range(GRIDT),
            lambda_GeneNet = lam_gn, kappa_GeneNet = kap_gn,
            lambda_beam = lam_b, kappa_beam = kap_b, alpha_beam = bf@alphaOpt,
            lambda_HS = lam_hs, kappa_HS = kappa_hat(spC, lam_hs, neff),
            K_eta = K_eta, eta0 = lf_rt$eta0,
            lambda_stab = st$pick[1], kappa_stab = kappa_hat(sp, st$pick[1], neff),
            omega_stab = omegas[match(st$pick[1], st$lambda)],
            lambda_stabp = st$pick[2], omega_stabp = omegas[match(st$pick[2], st$lambda)],
            kappa_x3 = kappa_hat(sp, 3 * lam, neff),
            kappa_x10 = kappa_hat(sp, 10 * lam, neff),
            kappa_x30 = kappa_hat(sp, 30 * lam, neff),
            m_liu = m_liu, m_true = m_true)

  out <- list(chk = chk, rank = rank, paths = paths, met3 = met3, tune = tune,
              re_info = re_info, pop_re = pop_re, stab = st$stab,
              cv_flat = c(range30 = diff(range(cv1)), range100 = diff(range(cv100))),
              ext = ext_new)
  if (r <= NOLD) {
    iv <- function(R) as.integer(round(1e5 * R[ut]))
    out$signed <- list(RT = iv(R_plain), vWP = iv(R_alt), vWPt = iv(R_altT),
                       GeneNet = iv(R_gn),
                       beam = as.integer(round(1e5 * btab[, "p_cor"])))
    if (par_ok && !is.null(pr)) out$signed$parcor_ridge <- iv(pr)
  }
  out$seconds <- proc.time()[3] - t0
  saveRDS(out, paste0(f_out, ".tmp"))
  file.rename(paste0(f_out, ".tmp"), f_out)
  "done"
}

## ---------------------------------------------------------------------------
## 3. The p = 400 block of sim_extensions.R, rebuilt
## ---------------------------------------------------------------------------
## Same seeds and the same order of random numbers as run_hd_task(): the
## graph and weights, then the data, then the CV folds. Nothing after the CV
## uses random numbers.
hd2_rep <- function(job) {
  f_out <- file.path(job$work, sprintf("hd2_p%d_%s_%d_rep%03d.rds", job$p,
                                       job$graph, job$n, job$rep))
  if (file.exists(f_out)) return("skip")
  old <- readRDS(file.path(job$work0, sprintf("hd_p%d_%s_%d_rep%03d.rds", job$p,
                                              job$graph, job$n, job$rep)))
  assign(".Random.seed", job$seed, envir = .GlobalEnv)
  p <- job$p; n <- job$n
  mod <- make_hd_model(job$graph, p)
  Sig <- mod$Sig; E <- mod$E; ut <- upper.tri(E); e <- E[ut]
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sig)
  X <- scale(X, center = TRUE, scale = FALSE)
  neff <- n - 1
  S <- crossprod(X) / neff
  sp <- spec_S(S)
  s_f <- .Random.seed
  lam <- cv_lambda(X, GRID30, type = "ridge")
  assign(".Random.seed", s_f, envir = .GlobalEnv)
  fold <- sample(rep_len(1:5, n))
  cv30 <- cv_curve(X, fold, GRID30); cv100 <- cv_curve(X, fold, GRID100)
  kap <- kappa_hat(sp, lam, neff)
  R_plain <- pcor_of(ridge_Q(sp, lam))
  R_kap <- pcor_of(solve(Sig + kap * diag(p)))
  P_true <- pcor_of(mod$Theta)
  RR <- richardson_Q(sp, neff, lam, cfac = 2)
  l2 <- lambda_for_kappa(sp, neff, 2 * kap, lower = lam)
  l2g <- grid_lambda2(sp, neff, lam, GRID30)
  RG <- richardson_Q2(sp, neff, lam, l2g)
  R_orc <- R_plain - (R_kap - P_true)
  R_orl <- R_plain; R_orl[!E] <- R_plain[!E] - R_kap[!E]; diag(R_orl) <- 1
  lam100 <- GRID100[which.min(cv100)]
  m_true <- sum(e)
  sc <- list(RT = abs(R_plain[ut]), RE = abs(pcor_of(RR$T)[ut]),
             RE_g30 = abs(pcor_of(RG$T)[ut]), ORC = abs(R_orc[ut]),
             ORC_leak = abs(R_orl[ut]),
             RT_cv100 = abs(pcor_of(ridge_Q(sp, lam100))[ut]))
  rank <- t(sapply(sc, rank_summary, e = e, m_true = m_true, m_liu = 0))
  k2 <- RR$kappa[2]
  out <- list(
    chk = c(lambda = lam - old$info[["lambda"]],
            orc_RT = rank["RT", "best"] - old$orc[["RT"]]),
    rank = rank,
    pop = rbind(RT = pop_summary(R_kap, E),
                RE = pop_summary(pcor_of(pop_extrap(Sig, kap, k2)), E)),
    re_info = c(lambda1 = lam, lambda2 = l2, ratio = l2 / lam,
                pd_guarantee = l2 <= 4 * lam,
                min_eig = min(eigen(RR$T, TRUE, only.values = TRUE)$values),
                lambda2_g30 = l2g, h_g30 = RG$h, kappa1 = kap, kappa2 = k2),
    cv30 = cv30, cv100 = cv100, lambda_cv100 = lam100,
    kappa_cv100 = kappa_hat(sp, lam100, neff))
  saveRDS(out, paste0(f_out, ".tmp"))
  file.rename(paste0(f_out, ".tmp"), f_out)
  "done"
}

## ---------------------------------------------------------------------------
## 4. Run
## ---------------------------------------------------------------------------
start_cluster <- function() {
  cl <- makePSOCKcluster(NCORE)
  clusterExport(cl, c("load_tools2", "CODE", "lib_extra", "p_sim", "n_grid",
                      "graphs", "GRID30", "GRID100", "GRIDT", "NOLD", "NPAR",
                      "stream_of", "cv_curve", "altridge_T", "daie_target",
                      "lfdr_fit", "vec_to_adj", "rank_summary", "path_summary",
                      "pop_summary", "richardson_Q2", "grid_lambda2",
                      "pop_extrap", "stab_select", "main2_rep", "hd2_rep"))
  clusterEvalQ(cl, { load_tools2(CODE, lib_extra); NULL })
  cl
}
tlog <- file.path(OUT, "results", "sim_main2_runtime.csv")
log_time <- function(part, t0, njobs) {
  row <- data.frame(part = part, jobs = njobs,
                    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    cores = NCORE, finished = format(Sys.time()))
  write.table(row, tlog, sep = ",", row.names = FALSE,
              col.names = !file.exists(tlog), append = file.exists(tlog))
}

if ("main" %in% args) {
  t0 <- Sys.time()
  kap <- read.csv(file.path(OUT, "results", "sim_kappa.csv"))
  jobs <- list()
  for (r in seq_len(NREP)) for (g in graphs) for (n in n_grid)
    jobs[[length(jobs) + 1]] <- list(graph = g, n = n, rep = r, work = WORK,
      work0 = WORK0, lam_fix = kap$lambda[kap$graph == g & kap$n == n])
  todo <- !file.exists(file.path(WORK, sprintf("main2_%s_%d_rep%03d.rds",
    sapply(jobs, `[[`, "graph"), sapply(jobs, `[[`, "n"), sapply(jobs, `[[`, "rep"))))
  jobs <- jobs[todo]
  cat(length(jobs), "main replicates to run\n")
  if (length(jobs)) {
    cl <- start_cluster()
    done <- parLapplyLB(cl, jobs, function(j)
      tryCatch(main2_rep(j), error = function(e)
        paste("ERROR", j$graph, j$n, j$rep, conditionMessage(e))))
    stopCluster(cl)
    err <- grep("^ERROR", unlist(done), value = TRUE)
    if (length(err)) writeLines(err, file.path(WORK, "errors_main.txt"))
    log_time("main design, 500 replicates", t0, length(jobs))
  }
  cat("main done\n")
}

if ("hd" %in% args) {
  t0 <- Sys.time()
  WORK_SAVE <- WORK
  jobs <- hd_jobs(400)                 # seeds of sim_extensions.R
  jobs <- lapply(jobs, function(j) { j$work <- WORK_SAVE; j$work0 <- WORK0; j })
  cl <- start_cluster()
  done <- parLapplyLB(cl, jobs, function(j)
    tryCatch(hd2_rep(j), error = function(e)
      paste("ERROR", j$graph, j$n, j$rep, conditionMessage(e))))
  stopCluster(cl)
  err <- grep("^ERROR", unlist(done), value = TRUE)
  if (length(err)) writeLines(err, file.path(WORK, "errors_hd.txt"))
  log_time("p = 400 block rebuilt", t0, length(jobs))
  cat("hd done\n")
}

writeLines(capture.output(sessionInfo()),
           file.path(OUT, "results", "sessionInfo_main2.txt"))
