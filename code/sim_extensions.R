###############################################################################
## Extensions of the simulation study for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Purpose: computations added to the study of sim_ridge_support.R.
##   main     reruns the main design (p = 100, 16 settings, 100 replicates)
##            with the same random numbers as sim_ridge_support.R and stores,
##            for every replicate, in addition to the published quantities:
##            the signed ridge partial correlations of every pair, the
##            population ridge network at the estimated effective penalty
##            (margin, leakage, best MCC), thresholded marginal correlations,
##            lFDR selection at cut-off 0.99, Ha and Sun (2014) and Liu (2013).
##   kgrid    reruns the effective-penalty grid of Table 1 and stores every
##            replicate (for Monte Carlo standard errors).
##   highdim  p >> n block: p = 400, n in {50, 100}, band, hub and random
##            graphs, a new graph with heterogeneous mixed-sign weights in
##            each replicate.
## The method code is read from sim_ridge_support.R (function definitions
## only), so the two scripts use the same code.
##
## Requirements: R >= 4.1; packages huge, glmnet, fdrtool, igraph, GGMridge,
##               SILGGM. An extra library folder can be given in the
##               environment variable RIDGE_RLIB.
## Usage  : Rscript sim_extensions.R [main] [kgrid] [highdim] [p1000]
##          (no argument: main, kgrid and highdim). Run from code/.
##          Every replicate is written to ../results/ext_work/ as soon as it
##          is finished; a restarted run skips the replicates already there.
##          sim_summaries.R turns the stored replicates into tables.
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("main", "kgrid", "highdim")
OUT   <- ".."
WORK  <- file.path(OUT, "results", "ext_work")
## number of worker processes (environment variable RIDGE_CORES, default 22)
NCORE <- min(as.integer(Sys.getenv("RIDGE_CORES", "22")), detectCores() - 1)
dir.create(WORK, showWarnings = FALSE, recursive = TRUE)
CODE  <- normalizePath(".")

## Load the packages and the method functions of sim_ridge_support.R. Only
## the function definitions are evaluated (the study itself is not run).
load_tools <- function(code_dir, lib_extra) {
  if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
  suppressPackageStartupMessages({
    library(huge); library(glmnet); library(fdrtool); library(igraph)
    library(GGMridge); library(SILGGM)
  })
  ex <- parse(file.path(code_dir, "sim_ridge_support.R"))
  for (e in ex) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function")))
      eval(e, envir = .GlobalEnv)
  }
  ## constants of the main script
  assign("MULT", c(0.3, 1, 3, 10), envir = .GlobalEnv)
  assign("NPERM", 99, envir = .GlobalEnv)
  assign("ALPHA", 0.1, envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG")
  invisible(TRUE)
}
load_tools(CODE, lib_extra)

p_sim  <- 100
n_grid <- c(50, 100, 200, 400)
graphs <- c("band", "random", "hub", "scale-free")

## ---------------------------------------------------------------------------
## 1. Helpers
## ---------------------------------------------------------------------------

## Run an expression that uses random numbers without changing the stream
## of the main replicate: the state is saved before and restored after.
keep_rng <- function(expr) {
  old <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", old, envir = .GlobalEnv))
  expr
}

## Population margin quantities of a partial-correlation matrix R:
## zeta = smallest |rho| over edges, ell = largest |rho| over non-edges.
pop_margin <- function(R, E) {
  ut <- upper.tri(E)
  s <- abs(R[ut]); e <- E[ut]
  c(zeta = min(s[e]), ell = max(s[!e]))
}

## Ha and Sun (2014): ridge inverse of the correlation matrix with the
## penalty chosen by 10-fold CV of p ridge regressions, Fisher z of the
## partial correlations, Efron's empirical null, Benjamini-Hochberg at the
## levels q. The settings follow the examples of the GGMridge package.
ha_sun_select <- function(X, q = c(0.1, 0.01)) {
  n <- nrow(X); p <- ncol(X); ut <- upper.tri(diag(p))
  sx <- scale(X, center = TRUE, scale = TRUE)
  lam_grid <- seq(0.1, 20, by = 0.1) * (n - 1)
  fit <- lambda.cv(x = sx, lambda = lam_grid, fold = 10L)
  lam <- fit$lambda[which.min(fit$spe)] / (n - 1)
  part <- -scaledMat(solve(lam * diag(p) + cor(X)))
  out <- lapply(q, function(qq) matrix(FALSE, p, p))
  ef <- tryCatch(getEfronp(z = transFisher(part[ut]), bins = 50L, maxQ = 13),
                 error = function(e) NULL)
  if (!is.null(ef)) {
    padj <- p.adjust(ef$correctp, method = "BH")
    for (k in seq_along(q)) {
      A <- matrix(FALSE, p, p); A[ut] <- padj < q[k]
      out[[k]] <- A | t(A)
    }
  }
  list(sel = out, lambda = lam, ok = !is.null(ef))
}

## Liu (2013): GGM estimation with FDR control, nodewise lasso version
## (method GFC_L of the SILGGM package), at FDR level alpha.
liu_select <- function(X, alpha = 0.1) {
  p <- ncol(X)
  fit <- tryCatch({
    capture.output(s <- SILGGM(X, method = "GFC_L", alpha = alpha))
    s }, error = function(e) NULL)
  if (is.null(fit)) return(list(sel = matrix(FALSE, p, p), ok = FALSE))
  A <- as.matrix(fit$global_decision[[1]]) != 0
  diag(A) <- FALSE
  list(sel = A | t(A), ok = TRUE)
}

## ---------------------------------------------------------------------------
## 2. Main design: one replicate with the extended output
## ---------------------------------------------------------------------------
## The first part repeats one_rep() of sim_ridge_support.R step by step, in
## the same order, so that the random numbers and the published quantities
## are identical. Everything added afterwards either uses no random numbers
## or runs inside keep_rng().
one_rep_ext <- function(n, Sig, Theta, E, D, lam_fix) {
  p <- ncol(Sig)
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sig)
  X <- scale(X, center = TRUE, scale = FALSE)
  neff <- n - 1
  S <- crossprod(X) / neff
  sp <- spec_S(S)
  ut <- upper.tri(E)
  grid <- exp(seq(log(1e-2), log(20), length.out = 30))

  lam <- cv_lambda(X, grid, type = "ridge")
  R_plain <- pcor_of(ridge_Q(sp, lam))
  kap <- kappa_hat(sp, lam, neff)
  RR <- richardson_Q(sp, neff, lam, cfac = 2)
  R_rich <- pcor_of(RR$T)
  lam_a <- cv_lambda(X, grid, type = "alt")
  R_alt <- pcor_of(altridge_Q(sp, lam_a))
  R_kap <- pcor_of(solve(Sig + kap * diag(p)))

  mps <- mps_select(X)
  gl <- glasso_ebic(X, S)
  hmb <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.01,
              verbose = FALSE)
  mb_sel <- as.matrix(huge.select(hmb, criterion = "ric",
                                  verbose = FALSE)$refit) != 0

  sel <- list(
    vWP_lFDR = lfdr_select(R_alt),
    RT_lFDR  = lfdr_select(R_plain),
    RE_lFDR  = lfdr_select(R_rich),
    RS_lasso = screen_refit(X, R_plain, dn = min(p - 1, ceiling(n / log(n)))),
    RT_MPS   = mps$A,
    GLASSO   = gl$sel,
    MB       = mb_sel)
  for (k in names(sel)) diag(sel[[k]]) <- FALSE
  met <- t(sapply(sel, sel_metrics, E = E))
  orc <- c(RT = oracle_mcc_threshold(R_plain, E),
           RE = oracle_mcc_threshold(R_rich, E),
           vWP = oracle_mcc_threshold(R_alt, E),
           GLASSO = oracle_mcc(gl$path, E),
           MB = oracle_mcc(lapply(hmb$path, function(A) as.matrix(A) != 0), E))
  sens <- rbind(mcc = sapply(mps$R, oracle_mcc_threshold, E = E),
                kappa = sapply(mps$lambda, kappa_hat, sp = sp, n = neff))
  rk <- rbind(RT = rank_metrics(R_plain, E), RE = rank_metrics(R_rich, E),
              vWP = rank_metrics(R_alt, E))
  dd <- factor(pmin(D[ut], 4), levels = 1:4)
  leak <- rbind(obs  = tapply(abs(R_plain[ut]), dd, mean),
                pred = tapply(abs(R_kap[ut]), dd, mean),
                rich = tapply(abs(R_rich[ut]), dd, mean))
  sig_min <- min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values)
  info <- c(lambda = lam, kappa = kap, kappa2 = RR$kappa[2],
            floor = lam + sig_min * max(0, p - neff + 1) / neff,
            maxdev = max(abs(R_plain - R_kap)[ut]))

  ## ---- additions without random numbers ----
  ## population ridge network at kappa_hat and at the nominal CV penalty
  R_lam <- pcor_of(solve(Sig + lam * diag(p)))
  m_lam <- pop_margin(R_lam, E)
  pop <- c(pop_margin(R_kap, E),
           mcc = oracle_mcc_threshold(R_kap, E),
           auc = unname(rank_metrics(R_kap, E)["AUC"]),
           zeta_lam = unname(m_lam["zeta"]), ell_lam = unname(m_lam["ell"]),
           mcc_lam = oracle_mcc_threshold(R_lam, E))
  ## the same on the penalty grid of the MPS selection
  pop_sens <- sapply(sens["kappa", ], function(k) {
    Rk <- pcor_of(solve(Sig + k * diag(p)))
    c(pop_margin(Rk, E), mcc = oracle_mcc_threshold(Rk, E))
  })
  ## marginal correlations: best attainable MCC and ranking summaries
  C <- cov2cor(S)
  marg <- c(mcc = oracle_mcc_threshold(C, E), rank_metrics(C, E))
  ## lFDR at cut-off 0.99
  sel2 <- list(vWP_lFDR99 = lfdr_select(R_alt, cut = 0.99),
               RT_lFDR99  = lfdr_select(R_plain, cut = 0.99))
  ## effective penalty at a fixed nominal penalty (same in all replicates)
  kap_fix <- kappa_hat(sp, lam_fix, neff)
  ## is the extrapolated matrix positive definite, with a positive diagonal?
  ev_rich <- min(eigen(RR$T, symmetric = TRUE, only.values = TRUE)$values)
  pd_rich <- c(min_eig = ev_rich, min_diag = min(diag(RR$T)))

  ## ---- additions with random numbers (stream left unchanged) ----
  hs <- keep_rng(ha_sun_select(X))
  liu <- keep_rng(liu_select(X, alpha = 0.1))
  sel2$HS_BH10 <- hs$sel[[1]]
  sel2$HS_BH01 <- hs$sel[[2]]
  sel2$Liu2013 <- liu$sel
  for (k in names(sel2)) diag(sel2[[k]]) <- FALSE
  met2 <- t(sapply(sel2, sel_metrics, E = E))

  list(met = met, rk = rk, leak = leak, orc = orc, sens = sens, info = info,
       pop = pop, pop_sens = pop_sens, marg = marg, met2 = met2,
       kap_fix = kap_fix, pd_rich = pd_rich,
       hs = c(lambda = hs$lambda, ok = hs$ok), liu_ok = liu$ok,
       ## signed pair-level values for the new leakage figure
       rho_hat = round(R_plain[ut], 6), rho_kap = round(R_kap[ut], 6))
}

## One chain = one random-number stream of one setting. As in run_streams()
## of the main script, stream i is nextRNGStream() applied i times to the
## setting seed and serves replicates i, i + 2, i + 4, ... in that order.
## Each replicate is saved together with the state of the stream after it,
## so an interrupted chain restarts where it stopped.
run_chain <- function(job) {
  gname <- job$graph; n <- job$n; st <- job$stream
  nrep <- 100; nstream <- 2
  mod <- make_model(gname, p_sim)
  set.seed(1000 + 10 * match(gname, graphs) + match(n, n_grid))
  s <- .Random.seed
  for (i in seq_len(st)) s <- parallel::nextRNGStream(s)
  assign(".Random.seed", s, envir = .GlobalEnv)
  for (r in seq(st, nrep, by = nstream)) {
    f <- file.path(job$work, sprintf("main_%s_%d_rep%03d.rds", gname, n, r))
    if (file.exists(f)) {
      assign(".Random.seed", readRDS(f)$seed_after, envir = .GlobalEnv)
      next
    }
    t0 <- proc.time()[3]
    res <- one_rep_ext(n, mod$Sig, mod$Theta, mod$E, mod$D, job$lam_fix)
    res$seconds <- proc.time()[3] - t0
    ## write to a temporary file first, so a stopped run leaves no broken file
    saveRDS(list(res = res, seed_after = .Random.seed), paste0(f, ".tmp"))
    file.rename(paste0(f, ".tmp"), f)
  }
  paste(gname, n, st)
}

## ---------------------------------------------------------------------------
## 3. Table 1 grid with every replicate stored
## ---------------------------------------------------------------------------
## Same code and seed as section 7 of sim_ridge_support.R; kh and dl of every
## replicate are kept, so that the Monte Carlo standard errors of the three
## columns and of the paired difference kappa_hat - kappa_n can be computed.
run_kgrid <- function() {
  set.seed(7)
  pk <- 200
  mk <- huge.generator(d = pk, graph = "band", g = 1, v = 0.4, u = 0.1,
                       verbose = FALSE)
  Sk <- as.matrix(mk$sigma); Ck <- chol(Sk)
  smin <- min(eigen(Sk, TRUE, only.values = TRUE)$values)
  kap_tab <- expand.grid(lambda = c(0.01, 0.1, 0.5, 1), n = c(50, 100, 200, 400, 800))
  rows <- list()
  for (i in seq_len(nrow(kap_tab))) {
    l <- kap_tab$lambda[i]; nn <- kap_tab$n[i]
    kh <- dl <- numeric(40)
    for (r in 1:40) {
      X <- matrix(rnorm(nn * pk), nn, pk) %*% Ck
      kh[r] <- kappa_hat(spec_S(crossprod(X) / nn), l, nn)
      Xm <- X[-1, , drop = FALSE]
      Qm <- solve(crossprod(Xm) / nn + l * diag(pk))
      dl[r] <- sum(Sk * Qm) / nn
    }
    rows[[i]] <- data.frame(lambda = l, n = nn, rep = 1:40, kappa_hat = kh,
                            kappa_mc_rep = l * (1 + dl),
                            floor = l + smin * max(0, pk - nn + 1) / nn)
  }
  do.call(rbind, rows)
}

## ---------------------------------------------------------------------------
## 4. p >> n block
## ---------------------------------------------------------------------------
## Graph from huge.generator (only its adjacency is used); edge weights with
## |theta_ij| ~ U[0.2, 0.5] and random signs; the diagonal is set as in
## huge, |smallest eigenvalue of the off-diagonal part| + 0.1 + u with
## u = 0.1; the covariance is rescaled to a correlation matrix. A new graph
## and new weights are drawn in every replicate.
make_hd_model <- function(graph, p) {
  A <- switch(graph,
    band   = huge.generator(d = p, graph = "band", g = 1, verbose = FALSE)$theta,
    random = huge.generator(d = p, graph = "random", prob = 3 / p,
                            verbose = FALSE)$theta,
    hub    = huge.generator(d = p, graph = "hub", g = p / 20,
                            verbose = FALSE)$theta)
  A <- as.matrix(A) != 0
  diag(A) <- FALSE
  W <- matrix(0, p, p)
  up <- which(upper.tri(A) & A)
  W[up] <- runif(length(up), 0.2, 0.5) * sample(c(-1, 1), length(up), replace = TRUE)
  W <- W + t(W)
  diag(W) <- abs(min(eigen(W, TRUE, only.values = TRUE)$values)) + 0.1 + 0.1
  Sig <- cov2cor(solve(W))
  Sig <- (Sig + t(Sig)) / 2
  list(Sig = Sig, Theta = solve(Sig), E = A)
}

## mps_select() of sim_ridge_support.R for p > n, written with the thin
## singular value decomposition of the n x p data matrix instead of the
## p x p eigendecomposition: with S = V diag(d^2) V' (rank r <= n),
## (S + l I)^{-1} = I / l + V diag(1 / (d^2 + l) - 1 / l) V'. The
## permutations use the random numbers in the same order, so the selection
## is the same as that of mps_select() up to rounding; it is only faster.
mps_select_thin <- function(X, mult = MULT, B = NPERM, alpha = ALPHA) {
  n <- nrow(X) - 1; p <- ncol(X); ut <- upper.tri(diag(p))
  ridge_pcors <- function(Y, lams) {
    s <- svd(Y / sqrt(n), nu = 0)
    keep <- s$d > 1e-10 * s$d[1]
    V <- s$v[, keep, drop = FALSE]; d2 <- s$d[keep]^2
    lapply(lams, function(l)
      pcor_of(diag(1 / l, p) + V %*% (t(V) * (1 / (d2 + l) - 1 / l))))
  }
  lams <- mult * sum(X^2) / n / p        # tr(S) / p, S = X'X / n
  Robs <- ridge_pcors(X, lams)
  Tobs <- sapply(Robs, function(R) max(abs(R[ut])))
  Tb <- matrix(0, B, length(lams))
  for (b in 1:B) {
    Xb <- scale(apply(X, 2, sample), center = TRUE, scale = FALSE)
    Tb[b, ] <- sapply(ridge_pcors(Xb, lams), function(R) max(abs(R[ut])))
  }
  k <- ceiling((1 - alpha / length(lams)) * (B + 1))
  A <- matrix(FALSE, p, p)
  tau <- numeric(length(lams))
  for (g in seq_along(lams)) {
    tau[g] <- sort(c(Tobs[g], Tb[, g]))[k]
    A <- A | (abs(Robs[[g]]) > tau[g])
  }
  diag(A) <- FALSE
  list(A = A, R = Robs, lambda = lams, tau = tau)
}

one_rep_hd <- function(graph, p, n, nmc = 5, full = TRUE) {
  mod <- make_hd_model(graph, p)
  Sig <- mod$Sig; E <- mod$E
  X <- matrix(rnorm(n * p), n, p) %*% chol(Sig)
  X <- scale(X, center = TRUE, scale = FALSE)
  neff <- n - 1
  S <- crossprod(X) / neff
  sp <- spec_S(S)
  ut <- upper.tri(E)
  grid <- exp(seq(log(1e-2), log(20), length.out = 30))
  t0 <- proc.time()[3]

  lam <- cv_lambda(X, grid, type = "ridge")
  R_plain <- pcor_of(ridge_Q(sp, lam))
  kap <- kappa_hat(sp, lam, neff)
  sig_ev <- eigen(Sig, symmetric = TRUE, only.values = TRUE)$values
  m0 <- max(0, p - neff + 1)
  floor_b <- lam + min(sig_ev) * m0 / neff
  ## sum of the m0 smallest eigenvalues (Ky Fan form of the floor)
  floor_kf <- lam + sum(sort(sig_ev)[seq_len(m0)]) / neff
  ## Monte Carlo kappa_n = lambda (1 + delta_n), delta_n = tr(Sigma E Q_-) / n,
  ## Q_- built from neff - 1 fresh N(0, Sigma) vectors with divisor neff
  Cs <- chol(Sig)
  dl <- sapply(seq_len(nmc), function(b) {
    Z <- matrix(rnorm((neff - 1) * p), neff - 1, p) %*% Cs
    sum(Sig * solve(crossprod(Z) / neff + lam * diag(p))) / neff
  })
  kap_mc <- lam * (1 + mean(dl))
  t_ridge <- proc.time()[3] - t0

  R_kap <- pcor_of(solve(Sig + kap * diag(p)))
  ## full = FALSE (used for p = 1000) keeps only the ridge parts; MPS, the
  ## graphical lasso and MB are too slow there for 100 replicates
  sel <- list(RT_lFDR = lfdr_select(R_plain))
  orc <- c(RT = oracle_mcc_threshold(R_plain, E), GLASSO = NA, MB = NA)
  t_mps <- t_gl <- t_mb <- NA
  if (full) {
    t1 <- proc.time()[3]
    mps <- if (p > n) mps_select_thin(X) else mps_select(X)
    t_mps <- proc.time()[3] - t1
    t1 <- proc.time()[3]
    gl <- glasso_ebic(X, S)
    t_gl <- proc.time()[3] - t1
    t1 <- proc.time()[3]
    hmb <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.01,
                verbose = FALSE)
    mb_sel <- as.matrix(huge.select(hmb, criterion = "ric",
                                    verbose = FALSE)$refit) != 0
    t_mb <- proc.time()[3] - t1
    sel <- c(sel, list(RT_MPS = mps$A, GLASSO = gl$sel, MB = mb_sel))
    orc["GLASSO"] <- oracle_mcc(gl$path, E)
    orc["MB"] <- oracle_mcc(lapply(hmb$path, function(A) as.matrix(A) != 0), E)
  }
  for (k in names(sel)) diag(sel[[k]]) <- FALSE
  met <- t(sapply(sel, sel_metrics, E = E))
  pc_true <- abs(pcor_of(mod$Theta)[ut][E[ut]])
  list(met = met,
       info = c(lambda = lam, kappa = kap, kappa_mc = kap_mc, floor = floor_b,
                floor_kyfan = floor_kf, n_edges = sum(E[ut]),
                pc_min = min(pc_true), pc_med = median(pc_true),
                pc_max = max(pc_true), cond = max(sig_ev) / min(sig_ev)),
       orc = orc,
       rk = rank_metrics(R_plain, E),
       pop = c(pop_margin(R_kap, E), mcc = oracle_mcc_threshold(R_kap, E)),
       seconds = c(ridge = t_ridge, mps = t_mps, glasso = t_gl, mb = t_mb))
}

## One task = one replicate; its stream is nextRNGStream() applied rep times
## to the setting seed, so every replicate is reproducible on its own.
run_hd_task <- function(job) {
  f <- file.path(job$work, sprintf("hd_p%d_%s_%d_rep%03d.rds", job$p,
                                   job$graph, job$n, job$rep))
  if (file.exists(f)) return("skip")
  assign(".Random.seed", job$seed, envir = .GlobalEnv)
  res <- tryCatch(one_rep_hd(job$graph, job$p, job$n, full = job$p < 1000),
                  error = function(e) list(error = conditionMessage(e)))
  saveRDS(res, paste0(f, ".tmp"))
  file.rename(paste0(f, ".tmp"), f)
  "done"
}

hd_jobs <- function(p, nrep = 100) {
  jobs <- list()
  hd_graphs <- c("band", "hub", "random")
  for (g in hd_graphs) for (n in c(50, 100)) {
    set.seed(3000 + 100 * (p == 1000) + 10 * match(g, hd_graphs) + n / 50)
    s <- .Random.seed
    for (r in seq_len(nrep)) {
      s <- parallel::nextRNGStream(s)
      jobs[[length(jobs) + 1]] <- list(graph = g, p = p, n = n, rep = r,
                                       seed = s, work = WORK)
    }
  }
  jobs
}

## ---------------------------------------------------------------------------
## 5. Run
## ---------------------------------------------------------------------------
start_cluster <- function() {
  cl <- makePSOCKcluster(NCORE)
  clusterExport(cl, c("load_tools", "CODE", "lib_extra", "p_sim", "n_grid",
                      "graphs", "keep_rng", "pop_margin", "ha_sun_select",
                      "liu_select", "one_rep_ext", "run_chain",
                      "make_hd_model", "mps_select_thin", "one_rep_hd", "run_hd_task",
                      "WORK"))
  clusterEvalQ(cl, { load_tools(CODE, lib_extra); Sys.setenv(OMP_NUM_THREADS = 1); NULL })
  cl
}
tlog <- file.path(OUT, "results", "sim_extensions_runtime.csv")
log_time <- function(part, t0) {
  row <- data.frame(part = part, minutes = as.numeric(difftime(Sys.time(), t0,
                                                               units = "mins")),
                    cores = NCORE, finished = format(Sys.time()))
  write.table(row, tlog, sep = ",", row.names = FALSE,
              col.names = !file.exists(tlog), append = file.exists(tlog))
}

if ("main" %in% args) {
  t0 <- Sys.time()
  kap <- read.csv(file.path(OUT, "results", "sim_kappa.csv"))
  jobs <- list()
  for (g in graphs) for (n in n_grid) for (st in 1:2)
    jobs[[length(jobs) + 1]] <- list(graph = g, n = n, stream = st, work = WORK,
      lam_fix = kap$lambda[kap$graph == g & kap$n == n])
  ## the large-n chains are the slowest; start them first
  jobs <- jobs[order(-sapply(jobs, `[[`, "n"))]
  cl <- start_cluster()
  done <- parLapplyLB(cl, jobs, run_chain)
  stopCluster(cl)
  log_time("main design rerun", t0)
  cat("main design done\n")
}

if ("kgrid" %in% args) {
  t0 <- Sys.time()
  kg <- run_kgrid()
  write.csv(kg, file.path(OUT, "results", "sim_kappa_grid_reps.csv"),
            row.names = FALSE)
  log_time("Table 1 grid rerun", t0)
  cat("kgrid done\n")
}

for (pp in c(400, 1000)) {
  if (!((pp == 400 && "highdim" %in% args) || (pp == 1000 && "p1000" %in% args)))
    next
  t0 <- Sys.time()
  jobs <- hd_jobs(pp)
  cl <- start_cluster()
  done <- parLapplyLB(cl, jobs, run_hd_task)
  stopCluster(cl)
  log_time(paste0("p >> n block, p = ", pp), t0)
  cat("high-dimensional block p =", pp, "done\n")
}

writeLines(capture.output(sessionInfo()),
           file.path(OUT, "results", "sessionInfo_extensions.txt"))
