###############################################################################
## Rescoring of stored replicates for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Purpose: four additions that reuse the random numbers of the stored runs.
##   designs  D1-D4 of sim_designs2.R: the best attainable MCC of the GLASSO
##            and MB paths is recomputed as in the main design, the largest
##            MCC over the graphs on each 30-point path. The paths are not
##            stored, so the data of every replicate are rebuilt from its
##            stream and only the two paths are refitted (in D3 also after
##            removing the top K principal components). The rebuilt data are
##            checked against the stored fingerprint, and the path-entry
##            score of sim_designs2.R, recomputed from the refitted paths,
##            must equal the stored value.
##            In D1 the ridge fits RT and vWP are also rebuilt (same CV folds)
##            for the tail-area rules below.
##   main     main design (sim_main2.R, 500 replicates): the RT and vWP fits
##            are rebuilt (same data, same CV folds) for the tail-area rules.
##   tail-area rules on the RT and vWP partial correlations:
##            BH_theo  Benjamini-Hochberg at 0.1 on Fisher-z p-values with a
##                     theoretical null N(0, 1/m0) centred at zero; m0 is the
##                     effective degrees of freedom of the estimator, 1 / the
##                     variance of atanh(rho_hat) under Sigma = I at the same
##                     n, p and penalty (Monte Carlo, part "null")
##            BH_n3    the same with m = n - 3 (the example of the paper's
##                     Proposition on theoretical nulls)
##            Fdr10    tail-area Fdr of fdrtool (empirical null) below 0.1
##            BH_emp   Benjamini-Hochberg at 0.1 on the fdrtool null p-values
##   constants  the crossover n of beta_n^(2) < beta_n, the chain factor of
##            the leakage bound, and Monte Carlo standard errors of the
##            margin columns of the leakage-and-noise table.
## Usage  : Rscript sim_rescore.R [null] [designs] [main] [summary] [constants]
##          (default: all, in this order), from code/. Every replicate is
##          written to ../results/ext_work4/ as soon as it is finished; a
##          restarted run skips the replicates already there. RIDGE_CORES
##          sets the number of worker processes (default 20).
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) args <- c("null", "designs", "main", "summary", "constants")
OUT   <- ".."
RES   <- normalizePath(file.path(OUT, "results"))
W0    <- file.path(RES, "ext_work")
W2    <- file.path(RES, "ext_work2")
W3    <- file.path(RES, "ext_work3")
W4    <- file.path(RES, "ext_work4")
dir.create(W4, showWarnings = FALSE)
CODE  <- normalizePath(".")
NCORE <- min(as.integer(Sys.getenv("RIDGE_CORES", "20")), detectCores() - 1)

## ---------------------------------------------------------------------------
## 1. Tools: function definitions of the earlier scripts, nothing else
## ---------------------------------------------------------------------------
load_tools4 <- function(code_dir, lib_extra) {
  Sys.setenv(OMP_NUM_THREADS = 1)
  if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
  suppressPackageStartupMessages({
    library(huge); library(fdrtool); library(igraph)
  })
  take <- function(file, only = NULL, values = NULL) {
    for (e in parse(file.path(code_dir, file))) {
      if (!(is.call(e) && identical(e[[1]], as.name("<-")))) next
      nm <- as.character(e[[2]])
      is_fun <- is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
      if ((is_fun && (is.null(only) || nm %in% only)) || nm %in% values)
        eval(e, envir = .GlobalEnv)
    }
  }
  take("sim_ridge_support.R")
  take("sim_designs2.R",
       c("make_d1_model", "make_d3_model", "make_main_model", "d4_raw",
         "normal_scores", "path_score", "rank_summ", "cached", "design_cells"),
       values = c("CV_GRID", "D4_A", "D4_B", "D4_W"))
  take("sim_main2.R", c("stream_of", "lfdr_fit", "vec_to_adj"),
       values = c("GRID30"))
  assign("MODEL_CACHE", new.env(), envir = .GlobalEnv)
  RNGkind("L'Ecuyer-CMRG")
  invisible(TRUE)
}
load_tools4(CODE, lib_extra)

graphs <- c("band", "random", "hub", "scale-free")
n_main <- c(50, 100, 200, 400)
P_MAIN <- 100
Q_LEVEL <- 0.1

## ---------------------------------------------------------------------------
## 2. Theoretical null of atanh(rho_hat) for a ridge-type estimator
## ---------------------------------------------------------------------------
## Under Sigma = I every population partial correlation is zero, at every
## penalty. The null variance of atanh(rho_hat) at (n, p, lambda) is the mean
## square over R0 null data sets and all pairs; m0 = 1 / variance is the
## effective degrees of freedom (for the unpenalized estimator with n > p + 1
## it is close to n - p - 1, for a large penalty it is far above n - 3,
## because the ridge partial correlations shrink towards zero).
null_table <- function(ns, grid, p, R0 = 50, seed) {
  set.seed(seed)
  ut <- upper.tri(diag(p))
  rows <- list()
  for (n in ns) {
    ss <- matrix(0, length(grid), 2)
    for (b in seq_len(R0)) {
      X <- scale(matrix(rnorm(n * p), n, p), center = TRUE, scale = FALSE)
      sp <- spec_S(crossprod(X) / (n - 1))
      for (g in seq_along(grid)) {
        z1 <- atanh(pcor_of(ridge_Q(sp, grid[g]))[ut])
        z2 <- atanh(pcor_of(altridge_Q(sp, grid[g]))[ut])
        ss[g, ] <- ss[g, ] + c(mean(z1^2), mean(z2^2)) / R0
      }
    }
    rows[[length(rows) + 1]] <- data.frame(n = n, lambda = grid,
      m0_RT = 1 / ss[, 1], m0_vWP = 1 / ss[, 2])
  }
  do.call(rbind, rows)
}
null_file <- file.path(W4, "null_m0.rds")
if ("null" %in% args && !file.exists(null_file)) {
  t0 <- Sys.time()
  nt <- list(main = null_table(n_main, GRID30, P_MAIN, seed = 20261101),
             d1 = null_table(c(50, 100, 200, 400, 800), CV_GRID, P_MAIN,
                             seed = 20261102))
  saveRDS(nt, null_file)
  cat("null table:", format(Sys.time() - t0), "\n")
}
NULLTAB <- if (file.exists(null_file)) readRDS(null_file) else NULL
m0_of <- function(tab, n, lambda, type) {
  i <- which(tab$n == n & abs(tab$lambda - lambda) < 1e-9 * lambda)
  if (length(i) != 1) stop("penalty not on the null grid")
  tab[[paste0("m0_", type)]][i]
}

## ---------------------------------------------------------------------------
## 3. Tail-area selection rules on one partial-correlation matrix
## ---------------------------------------------------------------------------
tail_rules <- function(R, n, m0, E, prefix) {
  p <- ncol(R); ut <- upper.tri(R)
  r <- R[ut]
  z <- atanh(pmin(abs(r), 1 - 1e-15))
  bh <- function(pv) p.adjust(pv, "BH") < Q_LEVEL
  lf <- lfdr_fit(R)
  sel <- list(BH_theo = bh(2 * pnorm(-z * sqrt(m0))),
              BH_n3   = bh(2 * pnorm(-z * sqrt(n - 3))),
              Fdr10   = lf$qval < Q_LEVEL,
              BH_emp  = bh(lf$pval),
              lFDR90  = 1 - lf$lfdr >= 0.9)
  met <- t(sapply(sel, function(v) sel_metrics(vec_to_adj(v, p), E)))
  rownames(met) <- paste(prefix, names(sel), sep = "_")
  met
}

## ---------------------------------------------------------------------------
## 4. Best attainable MCC along a path of graphs (main-design definition)
## ---------------------------------------------------------------------------
path_best <- function(path, E) {
  met <- t(sapply(path, function(A) { diag(A) <- FALSE; sel_metrics(A, E) }))
  j <- which.max(met[, "MCC"])
  c(best = unname(met[j, "MCC"]), tpr = unname(met[j, "TPR"]),
    fdp = unname(met[j, "FDP"]), size = unname(met[j, "NE"]),
    nsteps = length(path))
}
fit_paths <- function(X, E, C) {
  p <- ncol(X)
  hg <- huge(X, method = "glasso", nlambda = 30,
             lambda.min.ratio = if (p > 100) 0.05 else 0.01, verbose = FALSE)
  hm <- huge(X, method = "mb", nlambda = 30, lambda.min.ratio = 0.01,
             verbose = FALSE)
  pg <- lapply(hg$path, function(A) as.matrix(A) != 0)
  pm <- lapply(hm$path, function(A) as.matrix(A) != 0)
  ## the path-entry score of sim_designs2.R (for the check against storage)
  Omin <- as.matrix(hg$icov[[length(hg$icov)]])
  old <- c(GLASSO = unname(rank_summ(path_score(pg, abs(pcor_of(Omin))), E)["best"]),
           MB = unname(rank_summ(path_score(pm, abs(C)), E)["best"]))
  list(GLASSO = path_best(pg, E), MB = path_best(pm, E), old = old)
}

## ---------------------------------------------------------------------------
## 5. One replicate of D1-D4
## ---------------------------------------------------------------------------
## Same model construction and the same order of random numbers as one_rep2()
## of sim_designs2.R: the fixed models with their own seeds, then the stream
## of the replicate, then (D1) the model, then the data. In D1 the stream then
## continues as in analyse(): ridge CV folds, five leave-one-out draws of
## n - 2 rows, vWP CV folds.
designs_rep <- function(job) {
  f_out <- file.path(job$w4, paste0("dsg_", job$name, ".rds"))
  if (file.exists(f_out)) return("skip")
  stored <- readRDS(file.path(job$w2, paste0(job$name, ".rds")))
  if (!is.null(stored$error)) return("stored error")
  d <- job$design; cell <- job$cell
  if (d %in% c("d2", "d4"))
    mod <- cached(paste0("main_", cell$graph), make_main_model(cell$graph))
  if (d == "d3") mod <- cached(paste("d3", cell$graph, cell$p, cell$K),
                               make_d3_model(cell$graph, cell$p, cell$K))
  assign(".Random.seed", job$seed, envir = .GlobalEnv)
  n <- cell$n
  if (d == "d1") mod <- make_d1_model(cell$graph, 100)
  Cs <- chol(mod$Sig); p <- ncol(Cs)
  if (d == "d4") {
    X <- normal_scores(d4_raw(n, Cs, cell$dist))
  } else {
    X <- scale(matrix(rnorm(n * p), n, p) %*% Cs, center = TRUE, scale = FALSE)
  }
  E <- mod$E; neff <- n - 1
  S <- crossprod(X) / neff
  sp <- spec_S(S)
  chk <- c(fp = max(abs(X[1, 1:3] - stored$fp)))
  out <- list()
  if (d == "d1") {
    lam <- cv_lambda(X, CV_GRID, type = "ridge")
    for (b in 1:5) invisible(rnorm((neff - 1) * p))
    lam_v <- cv_lambda(X, CV_GRID, type = "alt")
    R_rt <- pcor_of(ridge_Q(sp, lam)); R_v <- pcor_of(altridge_Q(sp, lam_v))
    met <- rbind(
      tail_rules(R_rt, n, m0_of(job$nulltab, n, lam, "RT"), E, "RT"),
      tail_rules(R_v, n, m0_of(job$nulltab, n, lam_v, "vWP"), E, "vWP"),
      stored$a$met[c("RT_lFDR", "vWP0_lFDR", "MB", "Liu", "HS_BH10"), ])
    chk <- c(chk, lambda = lam - stored$a$info[["lambda"]],
             lam_vwp = lam_v - stored$a$info[["lam_vwp0"]],
             lfdr_RT = max(abs(met["RT_lFDR90", ] - stored$a$met["RT_lFDR", ])),
             lfdr_vWP = max(abs(met["vWP_lFDR90", ] - stored$a$met["vWP0_lFDR", ])))
    out$met <- met
    out$tune <- c(lambda = lam, lambda_vwp = lam_v,
                  m0_RT = m0_of(job$nulltab, n, lam, "RT"),
                  m0_vWP = m0_of(job$nulltab, n, lam_v, "vWP"))
  }
  pa <- fit_paths(X, E, cov2cor(S))
  out$paths <- rbind(GLASSO = pa$GLASSO, MB = pa$MB)
  out$paths_old_stored <- stored$a$rks[c("GLASSO", "MB"), "best"]
  chk <- c(chk, score_GL = unname(pa$old["GLASSO"] - stored$a$rks["GLASSO", "best"]),
           score_MB = unname(pa$old["MB"] - stored$a$rks["MB", "best"]))
  if (d == "d3") {
    V <- svd(X, nu = 0, nv = cell$K)$v
    Xr <- X - (X %*% V) %*% t(V)
    pr <- fit_paths(Xr, E, cov2cor(crossprod(Xr) / neff))
    out$paths_r <- rbind(GLASSO = pr$GLASSO, MB = pr$MB)
    out$paths_r_old_stored <- stored$r$rks[c("GLASSO", "MB"), "best"]
    chk <- c(chk, score_GL_r = unname(pr$old["GLASSO"] - stored$r$rks["GLASSO", "best"]),
             score_MB_r = unname(pr$old["MB"] - stored$r$rks["MB", "best"]))
  }
  out$chk <- chk
  saveRDS(out, paste0(f_out, ".tmp"))
  file.rename(paste0(f_out, ".tmp"), f_out)
  "done"
}

## jobs of D1-D4 with the seeds of make_jobs() in sim_designs2.R
designs_jobs <- function(nulltab) {
  jobs <- list()
  for (d in c("d1", "d2", "d3", "d4")) {
    cl <- design_cells(d)
    for (i in seq_len(nrow(cl))) {
      set.seed(cl$seed_base[i])
      s <- .Random.seed
      cell <- as.list(cl[i, setdiff(names(cl), c("name", "seed_base"))])
      r <- 0
      repeat {
        r <- r + 1
        s <- parallel::nextRNGStream(s)
        nm <- sprintf("%s_rep%03d", cl$name[i], r)
        if (!file.exists(file.path(W2, paste0(nm, ".rds")))) break
        jobs[[length(jobs) + 1]] <- list(design = d, cell = cell, seed = s,
          name = nm, w2 = W2, w4 = W4, nulltab = nulltab)
      }
    }
  }
  jobs
}

## ---------------------------------------------------------------------------
## 6. One replicate of the main design
## ---------------------------------------------------------------------------
## Same stream as main2_rep() of sim_main2.R: replicates 1-100 start from the
## states stored by sim_extensions.R, replicates 101-500 from new streams;
## then the data, the ridge CV folds and the vWP CV folds.
main_rep <- function(job) {
  g <- job$graph; n <- job$n; r <- job$rep
  f_out <- file.path(job$w4, sprintf("tail_main_%s_%d_rep%03d.rds", g, n, r))
  if (file.exists(f_out)) return("skip")
  gi <- match(g, graphs); ni <- match(n, n_main)
  mod <- make_model(g, P_MAIN)
  st <- readRDS(file.path(job$w3, sprintf("main2_%s_%d_rep%03d.rds", g, n, r)))
  if (r <= 100) {
    s0 <- if (r <= 2) stream_of(1000 + 10 * gi + ni, r)
          else readRDS(file.path(job$w0, sprintf("main_%s_%d_rep%03d.rds",
                                                 g, n, r - 2)))$seed_after
    ext <- readRDS(file.path(job$w0, sprintf("main_%s_%d_rep%03d.rds", g, n, r)))$res
  } else {
    s0 <- stream_of(5000 + 10 * gi + ni, r - 100)
    ext <- st$ext
  }
  assign(".Random.seed", s0, envir = .GlobalEnv)
  p <- P_MAIN
  X <- scale(matrix(rnorm(n * p), n, p) %*% chol(mod$Sig), center = TRUE,
             scale = FALSE)
  neff <- n - 1
  sp <- spec_S(crossprod(X) / neff)
  lam <- cv_lambda(X, GRID30, type = "ridge")
  lam_a <- cv_lambda(X, GRID30, type = "alt")
  R_rt <- pcor_of(ridge_Q(sp, lam)); R_v <- pcor_of(altridge_Q(sp, lam_a))
  E <- mod$E
  m0r <- m0_of(job$nulltab, n, lam, "RT"); m0v <- m0_of(job$nulltab, n, lam_a, "vWP")
  met <- rbind(tail_rules(R_rt, n, m0r, E, "RT"),
               tail_rules(R_v, n, m0v, E, "vWP"),
               ext$met[c("RT_lFDR", "vWP_lFDR", "MB"), ],
               Liu = ext$met2["Liu2013", ])
  chk <- c(lambda = lam - st$tune[["lambda"]],
           kappa = kappa_hat(sp, lam, neff) - st$tune[["kappa"]],
           lam_vwp = lam_a - st$tune[["lambda_vWP"]],
           Fdr_RT = max(abs(met["RT_Fdr10", ] - st$met3["RT_Fdr10", ])),
           BHemp_RT = max(abs(met["RT_BH_emp", ] - st$met3["RT_BH10", ])),
           Fdr_vWP = max(abs(met["vWP_Fdr10", ] - st$met3["vWP_Fdr10", ])),
           lfdr_RT = max(abs(met["RT_lFDR90", ] - ext$met["RT_lFDR", ])),
           lfdr_vWP = max(abs(met["vWP_lFDR90", ] - ext$met["vWP_lFDR", ])))
  out <- list(met = met, chk = chk,
              tune = c(lambda = lam, lambda_vwp = lam_a, m0_RT = m0r, m0_vWP = m0v))
  saveRDS(out, paste0(f_out, ".tmp"))
  file.rename(paste0(f_out, ".tmp"), f_out)
  "done"
}

## ---------------------------------------------------------------------------
## 7. Run the replicates
## ---------------------------------------------------------------------------
start_cluster4 <- function() {
  cl <- makePSOCKcluster(NCORE)
  clusterExport(cl, c("load_tools4", "CODE", "lib_extra", "graphs", "n_main",
                      "P_MAIN", "Q_LEVEL", "m0_of", "tail_rules", "path_best",
                      "fit_paths", "designs_rep", "main_rep"))
  clusterEvalQ(cl, { load_tools4(CODE, lib_extra); NULL })
  cl
}
run_jobs <- function(jobs, fun, label) {
  t0 <- Sys.time()
  cl <- start_cluster4()
  done <- parLapplyLB(cl, jobs, function(j)
    tryCatch(fun(j), error = function(e) paste("ERROR", conditionMessage(e))))
  stopCluster(cl)
  done <- unlist(done)
  err <- grep("^ERROR", done, value = TRUE)
  if (length(err)) writeLines(err, file.path(W4, paste0("errors_", label, ".txt")))
  cat(label, ":", sum(done == "done"), "done,", sum(done == "skip"), "skipped,",
      length(err), "errors,", format(Sys.time() - t0), "\n")
}

if ("designs" %in% args) {
  jobs <- designs_jobs(NULLTAB$d1)
  ## largest problems first
  sz <- sapply(jobs, function(j) j$cell$n * (if (is.null(j$cell$p)) 100 else j$cell$p))
  run_jobs(jobs[order(-sz)], designs_rep, "designs")
}
if ("main" %in% args) {
  jobs <- list()
  for (r in 1:500) for (g in graphs) for (n in n_main)
    jobs[[length(jobs) + 1]] <- list(graph = g, n = n, rep = r, w0 = W0, w3 = W3,
                                     w4 = W4, nulltab = NULLTAB$main)
  run_jobs(jobs, main_rep, "main")
}

## ---------------------------------------------------------------------------
## 8. Summaries
## ---------------------------------------------------------------------------
mse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
## Holm-adjusted paired z-tests of x - y over the settings of one comparison
holm_z <- function(dmean, dse) {
  pv <- 2 * pnorm(-abs(dmean / dse))
  p.adjust(pv, "holm")
}

cell_vars <- function(cell) {
  k <- strsplit(cell, "_")[[1]]
  switch(k[1],
    d1 = data.frame(design = "d1", graph = k[2], p = 100L, K = NA, dist = "gauss", n = as.integer(k[3])),
    d2 = data.frame(design = "d2", graph = k[2], p = 100L, K = NA, dist = "gauss", n = as.integer(k[3])),
    d3 = data.frame(design = "d3", graph = k[2], p = as.integer(k[3]), K = as.integer(k[4]),
                    dist = "gauss", n = as.integer(k[5])),
    d4 = data.frame(design = "d4", graph = k[2], p = 100L, K = NA, dist = k[3], n = as.integer(k[4])))
}

## selection summary of a list of per-replicate metric matrices, with the
## paired difference in MCC against a reference rule
sel_summary <- function(L, ref = "MB") {
  M <- simplify2array(L)                      # rule x stat x rep
  nr <- dim(M)[3]
  mcc <- M[, "MCC", , drop = FALSE][, 1, ]
  fw <- apply(M[, "FDP", , drop = FALSE] > 0, 1, mean)
  dd <- sweep(mcc, 2, mcc[ref, ])
  data.frame(rule = dimnames(M)[[1]], nrep = nr,
    MCC = rowMeans(mcc), seMCC = apply(mcc, 1, mse),
    TPR = rowMeans(M[, "TPR", ]), seTPR = apply(M[, "TPR", ], 1, mse),
    FDP = rowMeans(M[, "FDP", ]), seFDP = apply(M[, "FDP", ], 1, mse),
    FWER = fw, seFWER = sqrt(fw * (1 - fw) / nr), NE = rowMeans(M[, "NE", ]),
    dMCC_vs_MB = rowMeans(dd), se_dMCC_vs_MB = apply(dd, 1, mse),
    row.names = NULL)
}

if ("summary" %in% args) {
  ## ---- D1-D4 ----
  files <- list.files(W4, "^dsg_d[1-4]_.*_rep[0-9]+\\.rds$", full.names = TRUE)
  cells <- sub("^dsg_", "", sub("_rep[0-9]+\\.rds$", "", basename(files)))
  allchk <- list(); prow <- list(); frow <- list(); trow <- list()
  ## column names of the path summaries ("best", "tpr", ...), also for
  ## files written with suffixed names such as "best.MCC"
  fix_cols <- function(a) {
    for (k in intersect(c("paths", "paths_r"), names(a)))
      colnames(a[[k]]) <- sub("\\..*$", "", colnames(a[[k]]))
    a
  }
  read_dsg <- function(f) lapply(lapply(f, readRDS), fix_cols)
  for (cell in unique(cells)) {
    L <- read_dsg(files[cells == cell])
    cv <- cell_vars(cell)
    C <- do.call(rbind, lapply(L, `[[`, "chk"))
    allchk[[cell]] <- data.frame(cv, nrep = nrow(C), t(apply(abs(C), 2, max)))
    for (variant in if (cv$design == "d3") c("raw", "removed") else "raw") {
      pn <- if (variant == "raw") "paths" else "paths_r"
      po <- if (variant == "raw") "paths_old_stored" else "paths_r_old_stored"
      for (m in c("GLASSO", "MB")) {
        nw <- sapply(L, function(a) a[[pn]][m, "best"])
        od <- sapply(L, function(a) a[[po]][[m]])
        sz <- sapply(L, function(a) a[[pn]][m, "size"])
        prow[[length(prow) + 1]] <- data.frame(cv, variant = variant, method = m,
          nrep = length(nw), best_graph = mean(nw), se_graph = mse(nw),
          best_score_old = mean(od), se_score_old = mse(od),
          diff = mean(nw - od), se_diff = mse(nw - od),
          size_at_best = mean(sz), row.names = NULL)
      }
    }
    if (cv$design == "d1") {
      trow[[cell]] <- data.frame(cv, sel_summary(lapply(L, `[[`, "met")))
    }
  }
  ## the checks differ by design (D1 and D3 have more); fill with NA
  cols <- unique(unlist(lapply(allchk, names)))
  chkd <- do.call(rbind, lapply(allchk, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
  write.csv(chkd, file.path(RES, "sim_rescore_check_designs.csv"), row.names = FALSE)
  paths <- do.call(rbind, prow)
  write.csv(paths, file.path(RES, "sim_designs2_paths.csv"), row.names = FALSE)
  if (length(trow))
    write.csv(do.call(rbind, trow), file.path(RES, "sim_designs2_tailrules.csv"),
              row.names = FALSE)

  ## ---- D3: factor-adjusted comparison ----
  for (cell in unique(cells[grepl("^d3_", cells)])) {
    cv <- cell_vars(cell)
    S2 <- lapply(list.files(W2, paste0("^", cell, "_rep[0-9]+\\.rds$"),
                            full.names = TRUE), readRDS)
    L <- read_dsg(files[cells == cell])
    nr <- length(L)
    stopifnot(length(S2) == nr)
    vals <- list()
    for (variant in c("raw", "removed")) {
      part <- if (variant == "raw") "a" else "r"
      pn <- if (variant == "raw") "paths" else "paths_r"
      v <- list(
        RT_rank = sapply(S2, function(s) s[[part]]$rks["RT", "best"]),
        vWP_rank = sapply(S2, function(s) s[[part]]$rks["vWP0", "best"]),
        beam_rank = sapply(S2, function(s) s[[part]]$rks["beam", "best"]),
        GeneNet_rank = sapply(S2, function(s) s[[part]]$rks["GeneNet", "best"]),
        Marginal_rank = sapply(S2, function(s) s[[part]]$rks["Marginal", "best"]),
        GLASSO_path = sapply(L, function(a) a[[pn]]["GLASSO", "best"]),
        MB_path = sapply(L, function(a) a[[pn]]["MB", "best"]),
        GLASSO_EBIC_sel = sapply(S2, function(s) s[[part]]$met["GLASSO", "MCC"]),
        MB_RIC_sel = sapply(S2, function(s) s[[part]]$met["MB", "MCC"]),
        RT_lFDR_sel = sapply(S2, function(s) s[[part]]$met["RT_lFDR", "MCC"]),
        vWP_lFDR_sel = sapply(S2, function(s) s[[part]]$met["vWP0_lFDR", "MCC"]))
      fdp <- list(GLASSO_EBIC_sel = sapply(S2, function(s) s[[part]]$met["GLASSO", "FDP"]),
                  MB_RIC_sel = sapply(S2, function(s) s[[part]]$met["MB", "FDP"]),
                  RT_lFDR_sel = sapply(S2, function(s) s[[part]]$met["RT_lFDR", "FDP"]),
                  vWP_lFDR_sel = sapply(S2, function(s) s[[part]]$met["vWP0_lFDR", "FDP"]))
      tpr <- list(GLASSO_EBIC_sel = sapply(S2, function(s) s[[part]]$met["GLASSO", "TPR"]),
                  MB_RIC_sel = sapply(S2, function(s) s[[part]]$met["MB", "TPR"]),
                  RT_lFDR_sel = sapply(S2, function(s) s[[part]]$met["RT_lFDR", "TPR"]),
                  vWP_lFDR_sel = sapply(S2, function(s) s[[part]]$met["vWP0_lFDR", "TPR"]))
      vals[[variant]] <- v
      for (m in names(v)) frow[[length(frow) + 1]] <- data.frame(cv, nrep = nr,
        variant = variant, method = m, MCC = mean(v[[m]]), seMCC = mse(v[[m]]),
        TPR = if (m %in% names(tpr)) mean(tpr[[m]]) else NA,
        FDP = if (m %in% names(fdp)) mean(fdp[[m]]) else NA, row.names = NULL)
    }
    ## the fair comparisons, paired over replicates
    cmp <- list(
      "RT_rank(raw) - MB_path(removed)" = vals$raw$RT_rank - vals$removed$MB_path,
      "RT_rank(raw) - GLASSO_path(removed)" = vals$raw$RT_rank - vals$removed$GLASSO_path,
      "RT_rank(removed) - MB_path(removed)" = vals$removed$RT_rank - vals$removed$MB_path,
      "best ridge rank - MB_path(removed)" =
        pmax(vals$raw$RT_rank, vals$removed$RT_rank) - vals$removed$MB_path,
      "RT_lFDR(raw) - MB_RIC(removed)" = vals$raw$RT_lFDR_sel - vals$removed$MB_RIC_sel,
      "RT_lFDR(removed) - MB_RIC(removed)" = vals$removed$RT_lFDR_sel - vals$removed$MB_RIC_sel)
    for (m in names(cmp)) frow[[length(frow) + 1]] <- data.frame(cv, nrep = nr,
      variant = "paired difference", method = m, MCC = mean(cmp[[m]]),
      seMCC = mse(cmp[[m]]), TPR = NA, FDP = NA, row.names = NULL)
  }
  if (length(frow))
    write.csv(do.call(rbind, frow), file.path(RES, "sim_designs2_factor_fair.csv"),
              row.names = FALSE)

  ## ---- main design tail rules ----
  mrow <- list(); mchk <- list()
  for (g in graphs) for (n in n_main) {
    f <- file.path(W4, sprintf("tail_main_%s_%d_rep%03d.rds", g, n, 1:500))
    f <- f[file.exists(f)]
    if (length(f) < 2) next
    L <- lapply(f, readRDS)
    C <- do.call(rbind, lapply(L, `[[`, "chk"))
    mchk[[paste(g, n)]] <- data.frame(graph = g, n = n, nrep = length(L),
                                      t(apply(abs(C), 2, max)))
    TU <- do.call(rbind, lapply(L, `[[`, "tune"))
    mrow[[paste(g, n)]] <- data.frame(graph = g, n = n,
      sel_summary(lapply(L, `[[`, "met")),
      m0_RT_median = median(TU[, "m0_RT"]), m0_vWP_median = median(TU[, "m0_vWP"]))
  }
  if (length(mrow)) {
    write.csv(do.call(rbind, mchk), file.path(RES, "sim_rescore_check_main.csv"),
              row.names = FALSE)
    write.csv(do.call(rbind, mrow), file.path(RES, "sim_main2_tailrules.csv"),
              row.names = FALSE)
  }
  writeLines(capture.output(sessionInfo()),
             file.path(RES, "sessionInfo_rescore.txt"))
}

## ---------------------------------------------------------------------------
## 9. Constants
## ---------------------------------------------------------------------------
if ("constants" %in% args) {
  ## (a) beta_n against beta_n^(2), band model with p = 200 (Table 2), lambda = 1
  Sig200 <- make_model("band", 200)$Sig
  smax <- max(eigen(Sig200, symmetric = TRUE, only.values = TRUE)$values)
  beta1 <- function(n, p, s, l) {
    v <- 2 * p * s^2 / (n^2 * l^2) + p * s^3 / (n^3 * l^3)
    s / (n * l^2) + sqrt(3 * s / (4 * l^3)) * sqrt(v)
  }
  beta2 <- function(n, p, s, l) {
    w <- (s * (48 * p + 12 * p^2)^(1 / 4) / (n * l) + 2 * s^1.5 * sqrt(p) / (n^1.5 * l^1.5))^2
    3 * s / (n * l^2) + s^2 * sqrt(p) / (n^2 * l^3) + sqrt(3 * s / (4 * l^3)) * w
  }
  nn <- c(50, 100, 140:160, 200, 400, 800)
  cr <- data.frame(n = nn, p = 200, lambda = 1, sigma_max = smax,
                   beta_n = beta1(nn, 200, smax, 1), beta2_n = beta2(nn, 200, smax, 1))
  cr$beta2_smaller <- cr$beta2_n < cr$beta_n
  ## first integer n from which beta_n^(2) < beta_n holds for all larger n
  ## (checked up to 10^5)
  big <- 1:100000
  ok <- beta2(big, 200, smax, 1) < beta1(big, 200, smax, 1)
  n_star <- max(which(!ok)) + 1
  cr$n_crossover <- n_star
  root <- uniroot(function(n) beta2(n, 200, smax, 1) - beta1(n, 200, smax, 1),
                  c(100, 200), tol = 1e-10)$root
  cr$n_crossover_real <- root
  write.csv(cr, file.path(RES, "const_beta_crossover.csv"), row.names = FALSE)

  ## (b) chain: first bound of (13) over the actual entry at distance two,
  ## Theta = I + theta (adjacent pairs), p = 400, middle row
  ch <- list()
  for (th in c(0.4, 0.45)) {
    pc <- 400
    Th <- diag(pc); Th[abs(row(Th) - col(Th)) == 1] <- th
    et <- eigen(Th, symmetric = TRUE, only.values = TRUE)$values
    Sg <- solve(Th)
    for (k in c(0.1, 0.5, 1, 2, 3, 5)) {
      Tk <- solve(Sg + k * diag(pc)); i <- pc / 2
      vr <- (1 + k * max(et)) / (1 + k * min(et))
      q <- (sqrt(vr) - 1) / (sqrt(vr) + 1)
      b2 <- (1 + sqrt(vr))^2 / (2 * k * (1 + k * max(et))) * q^2
      dec <- abs(Tk[i, i + 3:8] / Tk[i, i + 2:7])
      ch[[length(ch) + 1]] <- data.frame(theta = th, p = pc, kappa = k,
        bound_d2 = b2, actual_d2 = abs(Tk[i, i + 2]), factor = b2 / abs(Tk[i, i + 2]),
        q_bar = q, observed_decay = mean(dec))
    }
  }
  write.csv(do.call(rbind, ch), file.path(RES, "const_chain_factor.csv"), row.names = FALSE)

  ## (c) MCSEs of the margin columns of the leakage-and-noise table
  sh <- function(x) sqrt(mean(x) * (1 - mean(x)) / length(x))
  tr <- list()
  for (g in graphs) for (n in n_main) {
    P1 <- sapply(1:100, function(r) readRDS(file.path(W0,
            sprintf("main_%s_%d_rep%03d.rds", g, n, r)))$res$pop[c("zeta", "ell")])
    P5 <- sapply(1:500, function(r) readRDS(file.path(W3,
            sprintf("main2_%s_%d_rep%03d.rds", g, n, r)))$pop_re["om0",
            c("zeta_RT", "ell_RT")])
    for (b in list(list(P1, 100), list(P5, 500))) {
      z <- b[[1]][1, ]; l <- b[[1]][2, ]
      tr[[length(tr) + 1]] <- data.frame(block = "main, p = 100", graph = g, n = n,
        reps = b[[2]], zeta = mean(z), se_zeta = mse(z), ell = mean(l), se_ell = mse(l),
        margin = mean(z - l), se_margin = mse(z - l),
        share_ell_gt_zeta = mean(l > z), se_share = sh(l > z))
    }
  }
  for (g in c("band", "hub", "random")) for (n in c(50, 100)) {
    f <- file.path(W0, sprintf("hd_p400_%s_%d_rep%03d.rds", g, n, 1:100))
    P <- sapply(f[file.exists(f)], function(x) readRDS(x)$pop[c("zeta", "ell")])
    z <- P[1, ]; l <- P[2, ]
    tr[[length(tr) + 1]] <- data.frame(block = "p = 400", graph = g, n = n,
      reps = ncol(P), zeta = mean(z), se_zeta = mse(z), ell = mean(l), se_ell = mse(l),
      margin = mean(z - l), se_margin = mse(z - l),
      share_ell_gt_zeta = mean(l > z), se_share = sh(l > z))
  }
  write.csv(do.call(rbind, tr), file.path(RES, "const_table3_mcse.csv"), row.names = FALSE)
  cat("constants written\n")
}
