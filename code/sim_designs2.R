###############################################################################
## Further simulation designs for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Purpose: four designs that are not in the main study.
##   d1   heterogeneous edge weights at p = 100: band, random, hub and
##        scale-free graphs, |theta_ij| ~ U[0.1, 0.6] with random signs, a
##        small diagonal offset (condition numbers about 20 to 160), a new
##        graph and new weights in every replicate; n in {50, ..., 800}.
##   d2   large n: the band and hub models of the main design with
##        n in {1600, 3200, 6400}.
##   d3   factor (spiked) covariance Sigma = B B' + Sigma_u with one or three
##        strong factors, the graph in Sigma_u^{-1}; p in {100, 200},
##        n in {100, 200, 400}; every method is also run after removing the
##        top K principal components.
##   d4   heavy tails and serial dependence: multivariate t with 5 degrees
##        of freedom and a GARCH(1,1) common-volatility design, analysed with
##        normal scores as in the data section; band and hub models.
##   kfit for d3 and d4 at fixed penalties: is the Monte Carlo mean of the
##        ridge estimator a ridge inverse c (Sigma + kappa I)^{-1}, and at
##        which kappa? Compared with the mean of kappa_hat and with the
##        leave-one-out kappa_n.
##   oracle (after the designs) splits the oracle bias correction of every
##        stored replicate into its leakage and edge-shrinkage parts.
## Every design compares the ridge rules with GLASSO, MB, Ha and Sun (2014),
## Liu (2013), GeneNet, beam and marginal correlations, both as selection
## rules and as rankings (best attainable MCC, and the share of true edges
## among the top |E| pairs).
##
## The method code of sim_ridge_support.R and sim_extensions.R is reused:
## only their function definitions are evaluated.
##
## Requirements: R >= 4.1; packages huge, glmnet, fdrtool, igraph, GGMridge,
##               SILGGM, GeneNet, beam. An extra library folder can be given
##               in the environment variable RIDGE_RLIB.
## Usage  : Rscript sim_designs2.R [d1] [d2] [d3] [d4] [kfit] [oracle] [quick]
##          (run_designs2.sh runs every part in the order used)
##          optional nrep=<k> (replicates per cell; default 100, 50 for d2)
##          and cells=<regex> (restrict to matching cell names).
##          Run from code/. Every replicate is written to
##          ../results/ext_work2/ as soon as it is finished; a restarted run
##          skips the replicates already there. Replicate r of a cell always
##          uses the same random-number stream, so raising nrep later only
##          adds replicates. sim_designs2_summary.R makes the tables.
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(key, default) {
  a <- grep(paste0("^", key, "="), args, value = TRUE)
  if (length(a)) sub(paste0("^", key, "="), "", a[1]) else default
}
QUICK <- "quick" %in% args
NREP_ARG <- as.integer(arg_val("nrep", NA))
CELLS <- arg_val("cells", ".")
OUT   <- ".."
WORK  <- arg_val("work", file.path(OUT, "results", "ext_work2"))
NCORE <- min(as.integer(Sys.getenv("RIDGE_CORES", "11")), detectCores() - 1)
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
CODE  <- normalizePath(".")

## Load packages and the function definitions of the two earlier scripts.
load_tools2 <- function(code_dir, lib_extra) {
  if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
  suppressPackageStartupMessages({
    library(huge); library(glmnet); library(fdrtool); library(igraph)
    library(GGMridge); library(SILGGM); library(GeneNet); library(beam)
  })
  take <- function(file, only = NULL) {
    for (e in parse(file.path(code_dir, file))) {
      if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
          identical(e[[3]][[1]], as.name("function")) &&
          (is.null(only) || as.character(e[[2]]) %in% only))
        eval(e, envir = .GlobalEnv)
    }
  }
  take("sim_ridge_support.R")
  take("sim_extensions.R", c("keep_rng", "pop_margin", "mps_select_thin"))
  assign("MULT", c(0.3, 1, 3, 10), envir = .GlobalEnv)
  assign("NPERM", 99, envir = .GlobalEnv)
  assign("ALPHA", 0.1, envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG")
  invisible(TRUE)
}
load_tools2(CODE, lib_extra)

## ---------------------------------------------------------------------------
## 1. Further estimators and selection rules
## ---------------------------------------------------------------------------

## Penalty grid of the likelihood cross-validation: 70 log-spaced values on
## [0.001, 50] (step ratio 1.17), finer and wider than in the main study.
CV_GRID <- exp(seq(log(1e-3), log(50), length.out = 70))

## van Wieringen-Peeters alternative ridge estimator with a scalar target
## T = tau I (ridgeP(S, lambda, type = "Alt", target = tau I) of
## rags2ridges): [ {lambda I + (S - lambda T)^2 / 4}^{1/2} + (S - lambda T) / 2
## ]^{-1}. With T = tau I it has the eigenvectors of S.
vwp_Q <- function(sp, lambda, tau) {
  s <- sp$val - lambda * tau
  sp$vec %*% (t(sp$vec) / (sqrt(lambda + s^2 / 4) + s / 2))
}

## Default target of rags2ridges, default.target(S) with type "DAIE": tau I,
## tau = mean of 1 / nu_j over the eigenvalues nu_j >= 1e-4 nu_max of S.
daie_tau <- function(sp) {
  v <- sp$val
  mean(1 / v[v >= 1e-4 * max(v)])
}

## K-fold cross-validated likelihood for the target version; as in
## rags2ridges the target is computed once, from the full data.
cv_vwp <- function(X, grid, tau, K = 5) {
  n <- nrow(X)
  fold <- sample(rep_len(1:K, n))
  loss <- matrix(0, K, length(grid))
  for (k in 1:K) {
    Xtr <- X[fold != k, , drop = FALSE]
    Xte <- X[fold == k, , drop = FALSE]
    sp <- spec_S(crossprod(Xtr) / nrow(Xtr))
    Ste <- crossprod(Xte) / nrow(Xte)
    ste_diag <- colSums(sp$vec * (Ste %*% sp$vec))
    for (g in seq_along(grid)) {
      s <- sp$val - grid[g] * tau
      d <- sqrt(grid[g] + s^2 / 4) + s / 2
      loss[k, g] <- sum(ste_diag / d) + sum(log(d))
    }
  }
  grid[which.min(colMeans(loss))]
}

## Ha and Sun (2014), as in sim_extensions.R, also returning the ridge
## partial correlations that it tests. For p = 200 or n >= 1600 the CV grid
## has step 0.5 instead of 0.1 (40 values on [0.1, 19.6]); with the full grid
## one fit takes about five minutes there.
ha_sun_fit <- function(X, q = 0.1) {
  n <- nrow(X); p <- ncol(X); ut <- upper.tri(diag(p))
  sx <- scale(X, center = TRUE, scale = TRUE)
  step <- if (p > 100 || n >= 1600) 0.5 else 0.1
  lam_grid <- seq(0.1, 20, by = step) * (n - 1)
  fit <- lambda.cv(x = sx, lambda = lam_grid, fold = 10L)
  lam <- fit$lambda[which.min(fit$spe)] / (n - 1)
  part <- -scaledMat(solve(lam * diag(p) + cor(X)))
  A <- matrix(FALSE, p, p)
  ef <- tryCatch(getEfronp(z = transFisher(part[ut]), bins = 50L, maxQ = 13),
                 error = function(e) NULL)
  if (!is.null(ef)) A[ut] <- p.adjust(ef$correctp, method = "BH") < q
  list(sel = A | t(A), R = part, lambda = lam, ok = !is.null(ef))
}

## Liu (2013), SILGGM method GFC_L at FDR level alpha; its test statistics
## give the ranking.
liu_fit <- function(X, alpha = 0.1) {
  p <- ncol(X)
  fit <- tryCatch({
    capture.output(s <- SILGGM(X, method = "GFC_L", alpha = alpha))
    s }, error = function(e) NULL)
  if (is.null(fit)) return(list(sel = matrix(FALSE, p, p),
                                Tst = matrix(0, p, p), ok = FALSE))
  A <- as.matrix(fit$global_decision[[1]]) != 0
  diag(A) <- FALSE
  Tst <- as.matrix(fit$T_stat); Tst[is.na(Tst)] <- 0
  list(sel = A | t(A), Tst = Tst, ok = TRUE)
}

## GeneNet (Schafer and Strimmer, 2005; Opgen-Rhein and Strimmer, 2007):
## shrinkage partial correlations (1 - l) R + l I with the analytic l of
## corpcor, edge tests of network.test.edges(). Selections: posterior edge
## probability 1 - lfdr >= 0.8 (the cut-off of extract.network()), >= 0.9
## (the cut-off of the lFDR rules), and Benjamini-Hochberg q-value < 0.1.
genenet_fit <- function(X) {
  p <- ncol(X)
  pc <- ggm.estimate.pcor(X, method = "static", verbose = FALSE)
  l <- attr(pc, "lambda")
  tt <- network.test.edges(pc, plot = FALSE, verbose = FALSE)
  mk <- function(keep) {
    A <- matrix(FALSE, p, p)
    A[cbind(tt$node1[keep], tt$node2[keep])] <- TRUE
    A | t(A)
  }
  R <- matrix(pc, p, p); diag(R) <- 1
  list(R = R, lambda = l,
       sel = list(GeneNet_p08 = mk(tt$prob >= 0.8), GeneNet_p09 = mk(tt$prob >= 0.9),
                  GeneNet_BH10 = mk(tt$qval < 0.1)))
}

## beam (Leday and Richardson, 2019), conditional dependencies with its
## default calibration, beam.select(thres = 0.1, method = "BH"). Its partial
## correlations are those of (R + {alpha / (1 - alpha)} I)^{-1}, R = cor(X),
## with the empirical-Bayes alpha, so the nominal ridge penalty on the
## correlation scale is alpha / (1 - alpha).
beam_fit <- function(X) {
  p <- ncol(X); ut <- upper.tri(diag(p))
  b <- beam(X, type = "conditional", return.only = c("cor", "prob"),
            verbose = FALSE)
  bs <- beam.select(b, thres = 0.1, method = "BH")
  R <- matrix(0, p, p); R[ut] <- b@table[, "p_cor"]; R <- R + t(R); diag(R) <- 1
  A <- matrix(FALSE, p, p)
  idx <- as.integer(rownames(bs@conditional))
  A[which(ut)[idx]] <- TRUE
  list(R = R, sel = A | t(A), alpha = b@alphaOpt)
}

## Thresholded marginal correlations: Benjamini-Hochberg at 0.1 on the
## t-tests of zero correlation.
marg_bh <- function(C, n, q = 0.1) {
  ut <- upper.tri(C)
  r <- C[ut]
  pv <- 2 * pt(-abs(r) * sqrt((n - 2) / (1 - r^2)), n - 2)
  A <- matrix(FALSE, nrow(C), ncol(C))
  A[ut] <- p.adjust(pv, "BH") < q
  A | t(A)
}

## Graphical lasso path tuned by EBIC (as glasso_ebic of the main script),
## also returning the score that ranks pairs along the path. For p = 200 the
## path ends at 0.05 times the largest penalty (as in the data analysis):
## with 0.01 one path takes about four minutes, and the path ending at 0.05
## already contains about 35 times as many edges as the graph.
glasso_fit <- function(X, nlambda = 30, gamma = 0.5) {
  n <- nrow(X); p <- ncol(X)
  S <- cor(X)
  h <- huge(X, method = "glasso", nlambda = nlambda,
            lambda.min.ratio = if (p > 100) 0.05 else 0.01, verbose = FALSE)
  eb <- sapply(seq_along(h$lambda), function(k) {
    Om <- as.matrix(h$icov[[k]])
    ne <- (sum(abs(Om) > 1e-10) - p) / 2
    ld <- as.numeric(determinant(Om, logarithm = TRUE)$modulus)
    n * (sum(S * Om) - ld) + ne * log(n) + 4 * gamma * ne * log(p)
  })
  path <- lapply(h$path, function(A) as.matrix(A) != 0)
  Omin <- as.matrix(h$icov[[length(h$icov)]])
  list(path = path, sel = path[[which.min(eb)]],
       score = path_score(path, abs(pcor_of(Omin))))
}

## Neighbourhood selection path, RIC tuning (as in the main script).
mb_fit <- function(X, C) {
  h <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.01,
            verbose = FALSE)
  sel <- as.matrix(huge.select(h, criterion = "ric", verbose = FALSE)$refit) != 0
  path <- lapply(h$path, function(A) as.matrix(A) != 0)
  list(path = path, sel = sel, score = path_score(path, abs(C)))
}

## A pair's score along a path: the number of path members that contain it
## (the paths are nested up to a few exceptions), ties broken by a small
## multiple of a secondary score.
path_score <- function(path, tie) {
  Reduce(`+`, lapply(path, function(A) A * 1)) + 1e-6 * tie / max(tie)
}

## Ranking summaries of a score matrix M (larger |M| = more likely an edge):
## best MCC over all thresholds, and the share of true edges among the top
## |E| pairs (TPR at matched size; this is also 1 - FDP).
rank_summ <- function(M, E) {
  ut <- upper.tri(E)
  s <- abs(M[ut]); e <- E[ut]
  o <- order(s, decreasing = TRUE)
  K <- sum(e)
  c(best = oracle_mcc_threshold(M, E), match = sum(e[o[seq_len(K)]]) / K)
}

## ---------------------------------------------------------------------------
## 2. The analysis of one data set
## ---------------------------------------------------------------------------
## X       centred n x p data
## E       true edge set
## Sigpop  covariance whose ridge inverses define the population ridge
##         network (Sigma, or the normal-scores covariance in d4)
## Ptarget partial correlations of the sparse precision matrix with support
##         E (used by the oracle bias correction)
## fresh   function(m) returning m new rows from the same law (for the
##         leave-one-out kappa_n), or NULL
## light   TRUE skips RT-MPS, Ha-Sun and Liu (used after removing principal
##         components)
analyse <- function(X, E, Sigpop, Ptarget, fresh = NULL, light = FALSE,
                    nmc = 5) {
  n <- nrow(X); p <- ncol(X); neff <- n - 1
  ut <- upper.tri(E)
  S <- crossprod(X) / neff
  sp <- spec_S(S)
  C <- cov2cor(S)
  tm <- c()
  t0 <- proc.time()[3]

  ## plain ridge at the likelihood cross-validated penalty
  lam <- cv_lambda(X, CV_GRID, type = "ridge")
  R_rt <- pcor_of(ridge_Q(sp, lam))
  kap <- kappa_hat(sp, lam, neff)
  ## leave-one-out kappa_n = lambda (1 + delta_n), delta_n = tr(Sigma E Q_-)
  ## / n, Q_- from neff - 1 new rows with divisor neff
  kap_mc <- NA
  if (!is.null(fresh)) {
    dl <- sapply(seq_len(nmc), function(b) {
      Z <- fresh(neff - 1)
      sum(Sigpop * solve(crossprod(Z) / neff + lam * diag(p))) / neff
    })
    kap_mc <- lam * (1 + mean(dl))
  }
  ## extrapolated ridge (kappa_2 = 2 kappa_1)
  RR <- richardson_Q(sp, neff, lam, cfac = 2)
  R_re <- pcor_of(RR$T)
  ## van Wieringen-Peeters, null target and rags2ridges default target
  lam_v0 <- cv_lambda(X, CV_GRID, type = "alt")
  R_v0 <- pcor_of(altridge_Q(sp, lam_v0))
  tau <- daie_tau(sp)
  lam_vd <- cv_vwp(X, CV_GRID, tau)
  R_vd <- pcor_of(vwp_Q(sp, lam_vd, tau))
  tm["ridge"] <- proc.time()[3] - t0

  sel <- list(RT_lFDR = lfdr_select(R_rt), RE_lFDR = lfdr_select(R_re),
              vWP0_lFDR = lfdr_select(R_v0), vWPd_lFDR = lfdr_select(R_vd))
  rk <- list(RT = R_rt, RE = R_re, vWP0 = R_v0, vWPd = R_vd, Marginal = C)
  sens <- NULL
  if (!light) {
    t0 <- proc.time()[3]
    mps <- if (p > neff) mps_select_thin(X) else mps_select(X)
    sel$RT_MPS <- mps$A
    kap_g <- sapply(mps$lambda, kappa_hat, sp = sp, n = neff)
    sens <- rbind(kappa = kap_g,
                  mcc = sapply(mps$R, oracle_mcc_threshold, E = E),
                  sapply(kap_g, function(k) {
                    Rk <- pcor_of(solve(Sigpop + k * diag(p)))
                    c(pop_margin(Rk, E), pop_mcc = oracle_mcc_threshold(Rk, E))
                  }))
    tm["mps"] <- proc.time()[3] - t0
  }
  t0 <- proc.time()[3]
  gl <- glasso_fit(X)
  tm["glasso"] <- proc.time()[3] - t0
  t0 <- proc.time()[3]
  mb <- mb_fit(X, C)
  tm["mb"] <- proc.time()[3] - t0
  sel$GLASSO <- gl$sel; sel$MB <- mb$sel
  rk$GLASSO <- gl$score; rk$MB <- mb$score
  t0 <- proc.time()[3]
  gn <- genenet_fit(X)
  bm <- beam_fit(X)
  tm["genenet_beam"] <- proc.time()[3] - t0
  sel <- c(sel, gn$sel, list(beam_BH10 = bm$sel, Marg_BH10 = marg_bh(C, n)))
  rk$GeneNet <- gn$R; rk$beam <- bm$R
  hs_ok <- liu_ok <- NA
  if (!light) {
    t0 <- proc.time()[3]
    hs <- ha_sun_fit(X)
    tm["hs"] <- proc.time()[3] - t0
    t0 <- proc.time()[3]
    liu <- liu_fit(X)
    tm["liu"] <- proc.time()[3] - t0
    sel$HS_BH10 <- hs$sel; sel$Liu <- liu$sel
    rk$HS <- hs$R; rk$Liu <- liu$Tst
    hs_ok <- hs$ok; liu_ok <- liu$ok
  }
  for (k in names(sel)) diag(sel[[k]]) <- FALSE
  met <- t(sapply(sel, sel_metrics, E = E))

  ## population ridge network at kappa_hat, and the oracle bias-corrected
  ## ranking rho_hat - (rho^(kappa_hat) - P(Theta)): sample noise without
  ## the population leakage and shrinkage
  R_kap <- pcor_of(solve(Sigpop + kap * diag(p)))
  rk$Oracle_corr <- R_rt - (R_kap - Ptarget)
  rks <- t(sapply(rk, rank_summ, E = E))
  ## (after removing principal components Sigpop is singular; then NA)
  P0 <- tryCatch(pcor_of(solve(Sigpop)),
                 error = function(e) matrix(NA, p, p))
  ## effective penalties of GeneNet and beam: ridge on the correlation
  ## matrix with nominal penalties l / (1 - l) and alpha / (1 - alpha)
  spC <- spec_S(C)
  kap_gn <- kappa_hat(spC, gn$lambda / (1 - gn$lambda), neff)
  kap_bm <- kappa_hat(spC, bm$alpha / (1 - bm$alpha), neff)
  Rgn <- pcor_of(solve(cov2cor(Sigpop) + kap_gn * diag(p)))
  ok0 <- !anyNA(P0)
  pop <- c(pop_margin(R_kap, E), pop_mcc = oracle_mcc_threshold(R_kap, E),
           pop_match = unname(rank_summ(R_kap, E)["match"]),
           zeta0 = if (ok0) unname(pop_margin(P0, E)["zeta"]) else NA,
           ell0 = if (ok0) unname(pop_margin(P0, E)["ell"]) else NA,
           pop0_mcc = if (ok0) oracle_mcc_threshold(P0, E) else NA,
           gn_zeta = unname(pop_margin(Rgn, E)["zeta"]),
           gn_ell = unname(pop_margin(Rgn, E)["ell"]),
           gn_pop_mcc = oracle_mcc_threshold(Rgn, E))
  info <- c(lambda = lam, kappa = kap, kappa_mc = kap_mc,
            kappa2 = RR$kappa[2], lam_vwp0 = lam_v0, lam_vwpd = lam_vd,
            tau = tau, gn_lambda = gn$lambda, kappa_gn = kap_gn,
            beam_alpha = bm$alpha, kappa_beam = kap_bm,
            hs_ok = hs_ok, liu_ok = liu_ok,
            maxdev = max(abs(R_rt - R_kap)[ut]))
  list(met = met, rks = rks, pop = pop, info = info, sens = sens, sec = tm)
}

## ---------------------------------------------------------------------------
## 3. Models
## ---------------------------------------------------------------------------

## d1: graph from huge.generator (adjacency only), |theta_ij| ~ U[0.1, 0.6]
## with random signs, diagonal |smallest eigenvalue of the off-diagonal
## part| + 0.1, rescaled to a correlation matrix. A new model per replicate.
make_d1_model <- function(graph, p) {
  A <- switch(graph,
    band   = huge.generator(d = p, graph = "band", g = 1, verbose = FALSE)$theta,
    random = huge.generator(d = p, graph = "random", prob = 3 / p,
                            verbose = FALSE)$theta,
    hub    = huge.generator(d = p, graph = "hub", g = 5, verbose = FALSE)$theta,
    "scale-free" = huge.generator(d = p, graph = "scale-free",
                                  verbose = FALSE)$theta)
  A <- as.matrix(A) != 0
  diag(A) <- FALSE
  W <- matrix(0, p, p)
  up <- which(upper.tri(A) & A)
  W[up] <- runif(length(up), 0.1, 0.6) * sample(c(-1, 1), length(up), replace = TRUE)
  W <- W + t(W)
  diag(W) <- abs(min(eigen(W, TRUE, only.values = TRUE)$values)) + 0.1
  Sig <- cov2cor(solve(W))
  Sig <- (Sig + t(Sig)) / 2
  list(Sig = Sig, E = A)
}

## d3: Sigma_u from the generator of the main design (band, bandwidth one;
## hub with p / 20 hubs), loadings of the first factor U[0.4, 0.8] (a common
## "market" factor), of the second and third N(0, 0.35^2); Sigma = B B' +
## Sigma_u, rescaled to a correlation matrix. One model per (graph, p, K).
make_d3_model <- function(graph, p, K) {
  set.seed(30000 + 100 * K + p + match(graph, c("band", "hub")))
  g <- switch(graph,
    band = huge.generator(d = p, graph = "band", g = 1, v = 0.4, u = 0.1,
                          verbose = FALSE),
    hub  = huge.generator(d = p, graph = "hub", g = p / 20, v = 0.3, u = 0.1,
                          verbose = FALSE))
  Su <- as.matrix(g$sigma)
  B <- cbind(runif(p, 0.4, 0.8), matrix(rnorm(p * 2, 0, 0.35), p, 2))[, 1:K,
                                                                     drop = FALSE]
  Sig <- cov2cor(tcrossprod(B) + Su)
  Sig <- (Sig + t(Sig)) / 2
  Thu <- solve(Su)
  E <- abs(Thu) > 1e-8; diag(E) <- FALSE
  ev <- eigen(Sig, TRUE)
  U <- ev$vectors[, 1:K, drop = FALSE]
  Pr <- diag(p) - tcrossprod(U)
  list(Sig = Sig, E = E, Ptarget = pcor_of(Thu), Sig_r = Pr %*% Sig %*% Pr,
       top_share = ev$values[1:K] / p, cond = max(ev$values) / min(ev$values))
}

## main-design models (band and hub of sim_ridge_support.R)
make_main_model <- function(graph) {
  m <- make_model(graph, 100)
  list(Sig = m$Sig, E = m$E, Ptarget = pcor_of(m$Theta))
}

## d4 data: X = sigma_t Z_t, Z_t ~ N(0, Sigma).
##   t5    sigma_t = (5 / W_t)^{1/2}, W_t ~ chi^2_5, independent rows
##   garch sigma_t^2 = omega + a eps_{t-1}^2 sigma_{t-1}^2 + b sigma_{t-1}^2,
##         eps_t = 1'Z_t / (1'Sigma 1)^{1/2} (the "market" shock), a = 0.15,
##         b = 0.80, omega = 0.05 (unit stationary variance), burn-in 500
D4_A <- 0.15; D4_B <- 0.80; D4_W <- 1 - D4_A - D4_B
d4_raw <- function(m, Cs, dist) {
  p <- ncol(Cs)
  Z <- matrix(rnorm(m * p), m, p) %*% Cs
  if (dist == "gauss") return(Z)
  if (dist == "t5") return(Z * sqrt(5 / rchisq(m, 5)))
  burn <- 500
  Zb <- matrix(rnorm(burn * p), burn, p) %*% Cs
  Zall <- rbind(Zb, Z)
  one <- rep(1, p)
  eps <- drop(Zall %*% one) / sqrt(sum(crossprod(Cs)))   # sum(Sigma) = 1'Sigma 1
  h <- numeric(burn + m); h[1] <- 1
  for (t in 2:(burn + m))
    h[t] <- D4_W + D4_A * eps[t - 1]^2 * h[t - 1] + D4_B * h[t - 1]
  (Zall * sqrt(h))[burn + seq_len(m), , drop = FALSE]
}
## within-sample normal scores Phi^{-1}{rank / (m + 1)}, centred
normal_scores <- function(X) {
  m <- nrow(X)
  Y <- apply(X, 2, function(x) qnorm(rank(x) / (m + 1)))
  scale(Y, center = TRUE, scale = FALSE)
}

## Population normal-scores covariance: Cov{Phi^{-1}(F_j(x_j))} with F_j
## the stationary marginal law, by Monte Carlo over 10^6 rows. For t5 the
## marginal is a t5 law; for garch F(x) = E Phi(x / sigma) over the
## stationary law of sigma_t, from a series of length 10^6.
d4_popcov <- function(Sig, dist, N = 1e6, chunk = 1e5) {
  p <- ncol(Sig); Cs <- chol(Sig)
  if (dist == "gauss") return(Sig)
  if (dist == "garch") {
    T0 <- 2e5 + 1000
    eps <- rnorm(T0); h <- numeric(T0); h[1] <- 1
    for (t in 2:T0) h[t] <- D4_W + D4_A * eps[t - 1]^2 * h[t - 1] + D4_B * h[t - 1]
    sig_s <- sqrt(h[-(1:1000)])
    xs <- seq(-25, 25, length.out = 5001)
    sub <- sig_s[seq(1, length(sig_s), length.out = 4000)]
    Fx <- sapply(xs, function(x) mean(pnorm(x / sub)))
    Fx <- pmin(pmax(Fx, 1e-12), 1 - 1e-12)
    g <- approxfun(xs, qnorm(Fx), rule = 2)
    draw_sig <- function(m) sample(sig_s, m, replace = TRUE)
  } else {
    g <- function(x) qnorm(pt(x, 5))
    draw_sig <- function(m) sqrt(5 / rchisq(m, 5))
  }
  M <- matrix(0, p, p); s1 <- numeric(p)
  for (k in seq_len(N / chunk)) {
    Z <- (matrix(rnorm(chunk * p), chunk, p) %*% Cs) * draw_sig(chunk)
    Y <- matrix(g(Z), chunk, p)
    M <- M + crossprod(Y); s1 <- s1 + colSums(Y)
  }
  V <- M / N - tcrossprod(s1 / N)
  (V + t(V)) / 2
}

## ---------------------------------------------------------------------------
## 4. Replicates
## ---------------------------------------------------------------------------

MODEL_CACHE <- new.env()
cached <- function(key, expr) {
  if (!exists(key, envir = MODEL_CACHE, inherits = FALSE))
    assign(key, expr, envir = MODEL_CACHE)
  get(key, envir = MODEL_CACHE)
}

one_rep2 <- function(job) {
  d <- job$design; cell <- job$cell
  ## fixed models are built (with their own seeds) before the stream is set
  if (d == "d2") mod <- cached(paste0("main_", cell$graph), make_main_model(cell$graph))
  if (d == "d3") mod <- cached(paste("d3", cell$graph, cell$p, cell$K),
                               make_d3_model(cell$graph, cell$p, cell$K))
  if (d == "d4") {
    mod <- cached(paste0("main_", cell$graph), make_main_model(cell$graph))
    Signs <- job$d4pop[[paste(cell$graph, cell$dist)]]
  }
  assign(".Random.seed", job$seed, envir = .GlobalEnv)
  n <- cell$n
  out <- list()
  if (d == "d1") {
    mod <- make_d1_model(cell$graph, 100)
    mod$Ptarget <- pcor_of(solve(mod$Sig))
    Cs <- chol(mod$Sig); p <- ncol(Cs)
    fresh <- function(m) matrix(rnorm(m * p), m, p) %*% Cs
    X <- scale(fresh(n), center = TRUE, scale = FALSE)
    ev <- eigen(mod$Sig, TRUE, only.values = TRUE)$values
    pc <- abs(mod$Ptarget[upper.tri(mod$E) & mod$E])
    out$model <- c(cond = max(ev) / min(ev), n_edges = length(pc),
                   pc_min = min(pc), pc_med = median(pc), pc_max = max(pc))
    out$fp <- X[1, 1:3]
    out$a <- analyse(X, mod$E, mod$Sig, mod$Ptarget, fresh)
  }
  if (d == "d2") {
    Cs <- chol(mod$Sig); p <- ncol(Cs)
    fresh <- function(m) matrix(rnorm(m * p), m, p) %*% Cs
    X <- scale(fresh(n), center = TRUE, scale = FALSE)
    out$fp <- X[1, 1:3]
    out$a <- analyse(X, mod$E, mod$Sig, mod$Ptarget, fresh)
  }
  if (d == "d3") {
    Cs <- chol(mod$Sig); p <- ncol(Cs)
    fresh <- function(m) matrix(rnorm(m * p), m, p) %*% Cs
    X <- scale(fresh(n), center = TRUE, scale = FALSE)
    out$fp <- X[1, 1:3]
    out$a <- analyse(X, mod$E, mod$Sig, mod$Ptarget, fresh)
    ## remove the top K sample principal components and analyse again; the
    ## population counterpart projects out the top K eigenvectors of Sigma
    V <- svd(X, nu = 0, nv = cell$K)$v
    Xr <- X - (X %*% V) %*% t(V)
    out$r <- analyse(Xr, mod$E, mod$Sig_r, mod$Ptarget, NULL, light = TRUE)
  }
  if (d == "d4") {
    Cs <- chol(mod$Sig); p <- ncol(Cs)
    fresh <- function(m) normal_scores(d4_raw(m, Cs, cell$dist))
    X <- fresh(n)
    out$fp <- X[1, 1:3]
    out$a <- analyse(X, mod$E, Signs, mod$Ptarget, fresh)
  }
  out$seed_end <- .Random.seed
  out
}

run_task2 <- function(job) {
  f <- file.path(job$work, paste0(job$name, ".rds"))
  if (file.exists(f)) return("skip")
  t0 <- proc.time()[3]
  res <- tryCatch(one_rep2(job), error = function(e) list(error = conditionMessage(e)))
  res$seconds <- proc.time()[3] - t0
  saveRDS(res, paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  "done"
}

## Cells and their seeds. Replicate r of a cell uses nextRNGStream() applied
## r times to the cell seed.
design_cells <- function(d) {
  cl <- switch(d,
    d1 = expand.grid(graph = c("band", "random", "hub", "scale-free"),
                     n = c(50, 100, 200, 400, 800), stringsAsFactors = FALSE),
    d2 = expand.grid(graph = c("band", "hub"), n = c(1600, 3200, 6400),
                     stringsAsFactors = FALSE),
    d3 = expand.grid(graph = c("band", "hub"), p = c(100, 200), K = c(1, 3),
                     n = c(100, 200, 400), stringsAsFactors = FALSE),
    d4 = expand.grid(graph = c("band", "hub"), dist = c("t5", "garch"),
                     n = c(100, 200, 400), stringsAsFactors = FALSE))
  cl$name <- apply(cl, 1, function(r) paste(c(d, gsub(" ", "", r)), collapse = "_"))
  cl$seed_base <- 20261031L + 1000L * as.integer(substr(d, 2, 2)) + seq_len(nrow(cl))
  cl
}
make_jobs <- function(d, nrep, d4pop = NULL) {
  cl <- design_cells(d)
  cl <- cl[grepl(CELLS, cl$name), , drop = FALSE]
  jobs <- list()
  for (i in seq_len(nrow(cl))) {
    set.seed(cl$seed_base[i])
    s <- .Random.seed
    cell <- as.list(cl[i, setdiff(names(cl), c("name", "seed_base"))])
    for (r in seq_len(nrep)) {
      s <- parallel::nextRNGStream(s)
      jobs[[length(jobs) + 1]] <- list(design = d, cell = cell, seed = s,
        name = sprintf("%s_rep%03d", cl$name[i], r), work = WORK, d4pop = d4pop)
    }
  }
  jobs
}

## ---------------------------------------------------------------------------
## 5. Is E Q a ridge inverse? (fixed penalties, d3 and d4 models)
## ---------------------------------------------------------------------------
## For a fixed lambda: the Monte Carlo mean of Q_lambda over R samples, the
## kappa that best fits it as c (Sigpop + kappa I)^{-1} (Frobenius norm,
## relative residual) and the kappa whose partial correlations best fit the
## mean of the sample partial correlations (root mean square over pairs);
## mean kappa_hat and the leave-one-out kappa_n. For d4 also on the raw
## (not normal-scored) data, with Sigpop = Sigma.
kfit_cell <- function(job) {
  f <- file.path(job$work, paste0(job$name, ".rds"))
  if (file.exists(f)) return("skip")
  cell <- job$cell
  if (cell$design == "d3") {
    mod <- make_d3_model(cell$graph, cell$p, cell$K)
    Cs <- chol(mod$Sig); Sigpop <- mod$Sig
    draw <- function(m) matrix(rnorm(m * ncol(Cs)), m, ncol(Cs)) %*% Cs
  } else {
    mod <- make_main_model(cell$graph)
    Cs <- chol(mod$Sig)
    if (cell$scores) {
      Sigpop <- job$d4pop[[paste(cell$graph, cell$dist)]]
      draw <- function(m) normal_scores(d4_raw(m, Cs, cell$dist))
    } else {
      Sigpop <- mod$Sig
      draw <- function(m) d4_raw(m, Cs, cell$dist)
    }
  }
  assign(".Random.seed", job$seed, envir = .GlobalEnv)
  p <- ncol(Cs); n <- cell$n; ut <- upper.tri(Sigpop)
  R <- job$R
  rows <- list()
  for (lam in job$lambdas) {
    Qm <- matrix(0, p, p); Pm <- numeric(sum(ut)); P2 <- Pm
    kh <- dl <- numeric(R)
    for (r in seq_len(R)) {
      X <- draw(n)
      X <- scale(X, center = TRUE, scale = FALSE)
      sp <- spec_S(crossprod(X) / (n - 1))
      Q <- ridge_Q(sp, lam)
      Qm <- Qm + Q / R
      pr <- pcor_of(Q)[ut]; Pm <- Pm + pr / R; P2 <- P2 + pr^2 / R
      kh[r] <- kappa_hat(sp, lam, n - 1)
      Z <- draw(n - 2)
      dl[r] <- sum(Sigpop * solve(crossprod(Z) / (n - 1) + lam * diag(p))) / (n - 1)
    }
    ev <- eigen(Sigpop, TRUE)
    fitQ <- function(k) {
      A <- ev$vectors %*% (t(ev$vectors) / (ev$values + k))
      cc <- sum(Qm * A) / sum(A * A)
      sqrt(sum((Qm - cc * A)^2) / sum(Qm^2))
    }
    fitP <- function(k) sqrt(mean((pcor_of(solve(Sigpop + k * diag(p)))[ut] - Pm)^2))
    oQ <- optimize(fitQ, c(1e-4, 100), tol = 1e-6)
    oP <- optimize(fitP, c(1e-4, 100), tol = 1e-6)
    kmc <- lam * (1 + mean(dl))
    mcse_pair <- sqrt(mean((P2 - Pm^2) / R))
    rows[[length(rows) + 1]] <- data.frame(job$cell, lambda = lam, R = R,
      kappa_hat = mean(kh), se_kappa_hat = sd(kh) / sqrt(R),
      kappa_loo = kmc, se_kappa_loo = lam * sd(dl) / sqrt(R),
      kappa_fitQ = oQ$minimum, relres_fitQ = oQ$objective,
      relres_at_kappa_hat = fitQ(mean(kh)),
      kappa_fitP = oP$minimum, rms_fitP = oP$objective,
      rms_at_kappa_hat = fitP(mean(kh)), rms_at_lambda = fitP(lam),
      mcse_pair_mean = mcse_pair)
  }
  res <- do.call(rbind, rows)
  saveRDS(res, paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  "done"
}

kfit_jobs <- function(d4pop, R) {
  cl <- rbind(
    data.frame(design = "d3", expand.grid(graph = c("band", "hub"), p = c(100, 200),
               K = c(1, 3), n = c(100, 200, 400), stringsAsFactors = FALSE),
               dist = "gauss", scores = FALSE, stringsAsFactors = FALSE),
    data.frame(design = "d4", p = 100, K = 0,
               expand.grid(graph = c("band", "hub"),
                           dist = c("gauss", "t5", "garch"), scores = c(TRUE, FALSE),
                           n = c(100, 200, 400), stringsAsFactors = FALSE),
               stringsAsFactors = FALSE))
  cl <- cl[!(cl$dist == "gauss" & !cl$scores & cl$design == "d4"), ]
  jobs <- list()
  for (i in seq_len(nrow(cl))) {
    nm <- paste("kfit", cl$design[i], cl$graph[i], cl$p[i], cl$K[i], cl$dist[i],
                if (cl$scores[i]) "ns" else "raw", cl$n[i], sep = "_")
    if (!grepl(CELLS, nm)) next
    set.seed(20269000L + i)
    jobs[[length(jobs) + 1]] <- list(cell = as.list(cl[i, ]), name = nm,
      seed = .Random.seed, R = R, lambdas = c(0.1, 0.5), work = WORK, d4pop = d4pop)
  }
  jobs
}

## ---------------------------------------------------------------------------
## 5b. Splitting the oracle bias correction (part "oracle")
## ---------------------------------------------------------------------------
## The oracle ranking rho_hat - (rho^(kappa_hat) - P(Theta)) removes both the
## leakage on non-edges and the shrinkage of the edges. This part rebuilds
## the data and the ridge fit of every stored replicate (same stream, same
## order of random numbers up to the cross-validation; checked against the
## stored penalty, kappa_hat and data fingerprint) and ranks
##   leak_only   rho_hat - rho^(kappa_hat) on non-edges (edges untouched)
##   shrink_only rho_hat + {P(Theta) - rho^(kappa_hat)} on edges only
##   full        the oracle ranking of analyse()
oracle_rep <- function(job) {
  f <- file.path(job$work, paste0("oracle_", job$name, ".rds"))
  if (file.exists(f)) return("skip")
  stored <- readRDS(file.path(job$work, paste0(job$name, ".rds")))
  if (!is.null(stored$error)) return("error")
  d <- job$design; cell <- job$cell
  if (d == "d2") mod <- cached(paste0("main_", cell$graph), make_main_model(cell$graph))
  if (d == "d3") mod <- cached(paste("d3", cell$graph, cell$p, cell$K),
                               make_d3_model(cell$graph, cell$p, cell$K))
  if (d == "d4") mod <- cached(paste0("main_", cell$graph), make_main_model(cell$graph))
  assign(".Random.seed", job$seed, envir = .GlobalEnv)
  n <- cell$n
  if (d == "d1") {
    mod <- make_d1_model(cell$graph, 100)
    mod$Ptarget <- pcor_of(solve(mod$Sig))
  }
  Cs <- chol(mod$Sig); p <- ncol(Cs)
  Sigpop <- mod$Sig
  if (d == "d4") {
    X <- normal_scores(d4_raw(n, Cs, cell$dist))
    Sigpop <- job$d4pop[[paste(cell$graph, cell$dist)]]
  } else {
    X <- scale(matrix(rnorm(n * p), n, p) %*% Cs, center = TRUE, scale = FALSE)
  }
  E <- mod$E; neff <- n - 1
  sp <- spec_S(crossprod(X) / neff)
  lam <- cv_lambda(X, CV_GRID, type = "ridge")
  kap <- kappa_hat(sp, lam, neff)
  R_rt <- pcor_of(ridge_Q(sp, lam))
  R_kap <- pcor_of(solve(Sigpop + kap * diag(p)))
  B <- R_kap - mod$Ptarget
  L <- R_kap; L[E] <- 0; diag(L) <- 0
  Sh <- -B; Sh[!E] <- 0
  rk <- list(RT = R_rt, leak_only = R_rt - L, shrink_only = R_rt + Sh,
             full = R_rt - B)
  chk <- c(lambda = lam - stored$a$info[["lambda"]],
           kappa = kap - stored$a$info[["kappa"]],
           fp = max(abs(X[1, 1:3] - stored$fp)))
  res <- list(rks = t(sapply(rk, rank_summ, E = E)), check = chk)
  saveRDS(res, paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  "done"
}

## ---------------------------------------------------------------------------
## 6. Run
## ---------------------------------------------------------------------------
start_cluster2 <- function() {
  cl <- makePSOCKcluster(NCORE)
  clusterExport(cl, c("load_tools2", "CODE", "lib_extra", "CV_GRID", "vwp_Q",
                      "daie_tau", "cv_vwp", "ha_sun_fit", "liu_fit",
                      "genenet_fit", "beam_fit", "marg_bh", "glasso_fit",
                      "mb_fit", "path_score", "rank_summ", "analyse",
                      "make_d1_model", "make_d3_model", "make_main_model",
                      "D4_A", "D4_B", "D4_W", "d4_raw", "normal_scores",
                      "MODEL_CACHE", "cached", "one_rep2", "run_task2",
                      "kfit_cell"))
  clusterEvalQ(cl, { load_tools2(CODE, lib_extra); Sys.setenv(OMP_NUM_THREADS = 1)
                     MODEL_CACHE <- new.env(); NULL })
  cl
}
tlog <- file.path(OUT, "results", "sim_designs2_runtime.csv")
log_time <- function(part, t0, njobs) {
  row <- data.frame(part = part, jobs = njobs,
                    minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")),
                    cores = NCORE, finished = format(Sys.time()))
  write.table(row, tlog, sep = ",", row.names = FALSE,
              col.names = !file.exists(tlog), append = file.exists(tlog))
}

## population normal-scores covariances of d4 (computed once, then read)
d4pop_file <- file.path(WORK, "d4_popcov.rds")
get_d4pop <- function() {
  if (file.exists(d4pop_file)) return(readRDS(d4pop_file))
  pc <- list()
  for (g in c("band", "hub")) {
    Sig <- make_main_model(g)$Sig
    for (dist in c("gauss", "t5", "garch")) {
      set.seed(20268000L + 10 * match(g, c("band", "hub")) +
               match(dist, c("gauss", "t5", "garch")))
      pc[[paste(g, dist)]] <- d4_popcov(Sig, dist)
    }
  }
  saveRDS(pc, d4pop_file)
  pc
}

designs <- intersect(c("d1", "d2", "d3", "d4"), args)
for (d in designs) {
  t0 <- Sys.time()
  nrep <- if (!is.na(NREP_ARG)) NREP_ARG else if (QUICK) 2 else if (d == "d2") 50 else 100
  jobs <- make_jobs(d, nrep, if (d == "d4") get_d4pop() else NULL)
  ## largest problems first
  sz <- sapply(jobs, function(j) j$cell$n * (if (is.null(j$cell$p)) 100 else j$cell$p))
  jobs <- jobs[order(-sz)]
  cl <- start_cluster2()
  done <- parLapplyLB(cl, jobs, run_task2)
  stopCluster(cl)
  log_time(d, t0, sum(unlist(done) == "done"))
  cat(d, "done:", sum(unlist(done) == "done"), "new replicates\n")
}
if ("oracle" %in% args) {
  t0 <- Sys.time()
  jobs <- list()
  for (d in c("d1", "d2", "d3", "d4")) {
    nrep <- if (!is.na(NREP_ARG)) NREP_ARG else if (d %in% c("d2")) 50 else 100
    jd <- make_jobs(d, nrep, if (d == "d4") get_d4pop() else NULL)
    have <- file.exists(file.path(WORK, paste0(sapply(jd, `[[`, "name"), ".rds")))
    jobs <- c(jobs, jd[have])
  }
  cl <- start_cluster2()
  clusterExport(cl, "oracle_rep")
  done <- parLapplyLB(cl, jobs, oracle_rep)
  stopCluster(cl)
  log_time("oracle split", t0, sum(unlist(done) == "done"))
  cat("oracle split done\n")
}
if ("kfit" %in% args) {
  t0 <- Sys.time()
  jobs <- kfit_jobs(get_d4pop(), R = if (QUICK) 5 else 200)
  cl <- start_cluster2()
  done <- parLapplyLB(cl, jobs, kfit_cell)
  stopCluster(cl)
  log_time("kfit", t0, sum(unlist(done) == "done"))
  cat("kfit done\n")
}

writeLines(capture.output(sessionInfo()),
           file.path(OUT, "results", "sessionInfo_designs2.txt"))
