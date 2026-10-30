###############################################################################
## Tables from the extended simulation runs (sim_extensions.R)
##
## Author : M. Arashi
## Purpose: reads the replicates stored in ../results/ext_work/ and writes
##   sim_extensions_check.csv  rerun against the published sim_results.rds
##   sim_margin.csv            population margin, leakage / noise split
##   sim_fwer.csv              share of replicates with a false edge
##   sim_competitors.csv       extra competitors on the main design
##   sim_mcse.csv              Monte Carlo standard errors, Tables 1-5
##   sim_claims.csv            paired differences behind the Section 5 text
##   sim_pairs_signed.csv      pair-level mean signed partial correlations
##   sim_leakage.csv           leakage summary by graph distance
##   sim_highdim.csv           p >> n block
## Usage  : Rscript sim_summaries.R   (run from code/)
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages({ library(huge); library(igraph) })
ex <- parse("sim_ridge_support.R")
for (e in ex) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) eval(e)
}
OUT  <- ".."
RES  <- file.path(OUT, "results")
WORK <- file.path(RES, "ext_work")
p_sim  <- 100
n_grid <- c(50, 100, 200, 400)
graphs <- c("band", "random", "hub", "scale-free")
MULT   <- c(0.3, 1, 3, 10)
keys   <- as.vector(t(outer(graphs, n_grid, paste, sep = "_")))
se     <- function(x) sd(x) / sqrt(length(x))
key_g  <- function(k) sub("_[0-9]+$", "", k)
key_n  <- function(k) as.integer(sub(".*_", "", k))

## ---------------------------------------------------------------------------
## Load the main-design replicates
## ---------------------------------------------------------------------------
reps <- lapply(keys, function(k) lapply(1:100, function(r)
  readRDS(file.path(WORK, sprintf("main_%s_rep%03d.rds", k, r)))$res))
names(reps) <- keys
old <- readRDS(file.path(RES, "sim_results.rds"))

## the rerun must reproduce the published replicates exactly
chk <- do.call(rbind, lapply(keys, function(k) {
  d <- sapply(c("met", "rk", "leak", "orc", "sens", "info"), function(q)
    max(sapply(1:100, function(r)
      max(abs(reps[[k]][[r]][[q]] - old[[k]][[r]][[q]]), na.rm = TRUE))))
  data.frame(setting = k, t(d))
}))
write.csv(chk, file.path(RES, "sim_extensions_check.csv"), row.names = FALSE)
cat("largest difference to sim_results.rds:", max(chk[, -1]), "\n")

## per-setting array helpers
arr_of <- function(k, q) simplify2array(lapply(reps[[k]], `[[`, q))
vec_of <- function(k, q) sapply(reps[[k]], `[[`, q)

## ---------------------------------------------------------------------------
## Task 1: population margin and the leakage / noise split
## ---------------------------------------------------------------------------
## the random and scale-free graphs depend on the generator kind: use the
## one of sim_ridge_support.R, otherwise set.seed() gives other graphs
RNGkind("L'Ecuyer-CMRG")
models <- lapply(graphs, make_model, p = p_sim); names(models) <- graphs
## check: the population ridge network of the rebuilt model must equal the
## stored one of the first replicate
for (k in keys) {
  m <- models[[key_g(k)]]
  rk1 <- pcor_of(solve(m$Sig + reps[[k]][[1]]$info["kappa"] * diag(p_sim)))
  stopifnot(max(abs(rk1[upper.tri(rk1)] - reps[[k]][[1]]$rho_kap)) < 1e-5)
}
pop_at <- function(Sig, E, kappa) {
  R <- pcor_of(solve(Sig + kappa * diag(ncol(Sig))))
  ut <- upper.tri(E); s <- abs(R[ut]); e <- E[ut]
  c(zeta = min(s[e]), ell = max(s[!e]))
}
marg <- do.call(rbind, lapply(keys, function(k) {
  g <- key_g(k); m <- models[[g]]
  I <- vec_of(k, "info"); P <- vec_of(k, "pop"); O <- vec_of(k, "orc")
  RK <- arr_of(k, "rk")
  at_k <- pop_at(m$Sig, m$E, mean(I["kappa", ]))
  at_l <- pop_at(m$Sig, m$E, mean(I["lambda", ]))
  base <- data.frame(graph = g, n = key_n(k))
  rows <- data.frame(base, penalty = c("kappa_hat at CV", "lambda_CV (nominal)"),
    kappa = c(mean(I["kappa", ]), mean(I["lambda", ])),
    zeta = c(at_k["zeta"], at_l["zeta"]), ell = c(at_k["ell"], at_l["ell"]),
    margin = c(at_k["zeta"] - at_k["ell"], at_l["zeta"] - at_l["ell"]),
    share_rep_margin_pos = c(mean(P["zeta", ] > P["ell", ]),
                             mean(P["zeta_lam", ] > P["ell_lam", ])),
    mcc_pop = c(mean(P["mcc", ]), mean(P["mcc_lam", ])),
    se_mcc_pop = c(se(P["mcc", ]), se(P["mcc_lam", ])),
    mcc_sample = c(mean(O["RT", ]), NA), se_mcc_sample = c(se(O["RT", ]), NA),
    sep_pop = c(mean(P["zeta", ] > P["ell", ]), NA),
    sep_sample = c(mean(RK["RT", "MARGIN", ] > 0), NA))
  ## the four penalties of the MPS grid (best MCC of the sample ranking from
  ## the published sens matrix, of the population ranking from pop_sens)
  SE <- arr_of(k, "sens"); PS <- arr_of(k, "pop_sens")
  grid_rows <- data.frame(base, penalty = paste0("grid m = ", MULT),
    kappa = rowMeans(SE["kappa", , ]),
    zeta = rowMeans(PS["zeta", , ]), ell = rowMeans(PS["ell", , ]),
    margin = rowMeans(PS["zeta", , ] - PS["ell", , ]),
    share_rep_margin_pos = rowMeans(PS["zeta", , ] > PS["ell", , ]),
    mcc_pop = rowMeans(PS["mcc", , ]), se_mcc_pop = apply(PS["mcc", , ], 1, se),
    mcc_sample = rowMeans(SE["mcc", , ]), se_mcc_sample = apply(SE["mcc", , ], 1, se),
    sep_pop = rowMeans(PS["zeta", , ] > PS["ell", , ]), sep_sample = NA)
  rbind(rows, grid_rows)
}))
## leakage error = 1 - best MCC of the population ranking; noise error = the
## further loss of the sample ranking
marg$leak_error  <- 1 - marg$mcc_pop
marg$noise_error <- marg$mcc_pop - marg$mcc_sample
write.csv(marg, file.path(RES, "sim_margin.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Task 2 and 6: error rates and extra competitors
## ---------------------------------------------------------------------------
met_rows <- function(k, q, design = "main") {
  A <- arr_of(k, q)
  fp <- round(A[, "FDP", ] * A[, "NE", ])
  data.frame(design = design, p = p_sim, graph = key_g(k), n = key_n(k),
             method = dimnames(A)[[1]],
             TPR = apply(A[, "TPR", ], 1, mean), FDP = apply(A[, "FDP", ], 1, mean),
             MCC = apply(A[, "MCC", ], 1, mean), seMCC = apply(A[, "MCC", ], 1, se),
             NE = apply(A[, "NE", ], 1, mean),
             FWER = rowMeans(fp > 0), seFWER = sqrt(rowMeans(fp > 0) *
                                                    (1 - rowMeans(fp > 0)) / 100),
             mean_FP = rowMeans(fp), row.names = NULL)
}
fwer <- do.call(rbind, c(lapply(keys, met_rows, q = "met"),
                         lapply(keys, met_rows, q = "met2")))
comp <- do.call(rbind, lapply(keys, met_rows, q = "met2"))
## thresholded marginal correlations: best attainable MCC (paired against RT)
mc <- do.call(rbind, lapply(keys, function(k) {
  M <- vec_of(k, "marg"); O <- vec_of(k, "orc")
  d <- M["mcc", ] - O["RT", ]
  data.frame(design = "main", p = p_sim, graph = key_g(k), n = key_n(k),
             method = "Marginal correlation (best threshold)",
             best_MCC = mean(M["mcc", ]), se_best_MCC = se(M["mcc", ]),
             AUC = mean(M["AUC", ]), sep = mean(M["MARGIN", ] > 0),
             RT_best_MCC = mean(O["RT", ]), diff_vs_RT = mean(d), se_diff = se(d))
}))
hs_fail <- sapply(keys, function(k) sum(vec_of(k, "hs")["ok", ] == 0))
liu_fail <- sapply(keys, function(k) sum(!vec_of(k, "liu_ok")))
comp$failures <- NA
comp$failures[grepl("HS", comp$method)] <- rep(hs_fail, each = 2)
comp$failures[comp$method == "Liu2013"] <- liu_fail
comp_all <- merge(comp, mc, all = TRUE, sort = FALSE)
write.csv(comp_all, file.path(RES, "sim_competitors.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Task 3: Monte Carlo standard errors of Tables 1-5
## ---------------------------------------------------------------------------
mc_rows <- list()
add <- function(table, graph, n, column, value, mcse, extra = NA)
  mc_rows[[length(mc_rows) + 1]] <<- data.frame(table = table, graph = graph,
    n = n, column = column, value = value, mcse = mcse, note = extra)

## Table 1 (rerun with every replicate kept)
kg <- read.csv(file.path(RES, "sim_kappa_grid_reps.csv"))
kg_pub <- read.csv(file.path(RES, "sim_kappa_grid.csv"))
for (i in seq_len(nrow(kg_pub))) {
  l <- kg_pub$lambda[i]; nn <- kg_pub$n[i]
  d <- kg[kg$lambda == l & kg$n == nn, ]
  stopifnot(abs(mean(d$kappa_hat) - kg_pub$kappa_hat[i]) < 1e-12,
            abs(mean(d$kappa_mc_rep) - kg_pub$kappa_mc[i]) < 1e-12)
  lab <- paste0("lambda=", l)
  add("Table 1", "band p=200", nn, paste(lab, "MC kappa_n"), mean(d$kappa_mc_rep), se(d$kappa_mc_rep))
  add("Table 1", "band p=200", nn, paste(lab, "est."), mean(d$kappa_hat), se(d$kappa_hat))
  add("Table 1", "band p=200", nn, paste(lab, "bound"), d$floor[1], 0, "deterministic")
  dd <- d$kappa_hat - d$kappa_mc_rep
  add("Table 1", "band p=200", nn, paste(lab, "est. - MC (paired)"), mean(dd), se(dd))
  add("Table 1", "band p=200", nn, paste(lab, "(est. - MC) / MC"),
      mean(dd) / mean(d$kappa_mc_rep), se(dd) / mean(d$kappa_mc_rep))
}
for (k in keys) {
  g <- key_g(k); n <- key_n(k)
  I <- vec_of(k, "info"); kf <- vec_of(k, "kap_fix")
  add("Table 2", g, n, "lambda_CV", mean(I["lambda", ]), se(I["lambda", ]))
  add("Table 2", g, n, "kappa_hat", mean(I["kappa", ]), se(I["kappa", ]))
  add("Table 2", g, n, "sd(kappa_hat)", sd(I["kappa", ]), sd(I["kappa", ]) / sqrt(2 * 99),
      "mcse of an sd, normal approximation")
  add("Table 2", g, n, "sd(kappa_hat) at fixed lambda", sd(kf), sd(kf) / sqrt(2 * 99),
      "lambda fixed at the mean lambda_CV of the setting")
  add("Table 2", g, n, "distinct lambda_CV values", length(unique(round(I["lambda", ], 10))), NA)
  add("Table 2", g, n, "bound", mean(I["floor", ]), se(I["floor", ]))
  add("Table 2", g, n, "max dev.", mean(I["maxdev", ]), se(I["maxdev", ]))
  add("Table 2", g, n, "kappa_hat / lambda_CV", mean(I["kappa", ] / I["lambda", ]),
      se(I["kappa", ] / I["lambda", ]))
  A <- arr_of(k, "met")
  for (m in dimnames(A)[[1]]) {
    add("Table 3", g, n, paste(m, "MCC"), mean(A[m, "MCC", ]), se(A[m, "MCC", ]))
    add("Table 4", g, n, paste(m, "TPR"), mean(A[m, "TPR", ]), se(A[m, "TPR", ]))
    add("Table 4", g, n, paste(m, "FDP"), mean(A[m, "FDP", ]), se(A[m, "FDP", ]))
    add("Table 4", g, n, paste(m, "NE"), mean(A[m, "NE", ]), se(A[m, "NE", ]))
  }
  O <- vec_of(k, "orc"); RK <- arr_of(k, "rk")
  for (m in rownames(O)) add("Table 5", g, n, paste(m, "best MCC"), mean(O[m, ]), se(O[m, ]))
  for (m in dimnames(RK)[[1]]) add("Table 5", g, n, paste(m, "AUC"), mean(RK[m, "AUC", ]),
                                   se(RK[m, "AUC", ]))
  ps <- mean(RK["RT", "MARGIN", ] > 0)
  add("Table 5", g, n, "Sep. (share of exact separation, RT)", ps, sqrt(ps * (1 - ps) / 100),
      "binomial")
}
mcse <- do.call(rbind, mc_rows)
write.csv(mcse, file.path(RES, "sim_mcse.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Task 3: the comparisons behind the sentences of Section 5
## ---------------------------------------------------------------------------
cl_rows <- list()
claim <- function(id, text, setting, quantity, d, value = mean(d), s = se(d))
  cl_rows[[length(cl_rows) + 1]] <<- data.frame(claim = id, sentence = text,
    setting = setting, quantity = quantity, estimate = value, se = s,
    z = ifelse(is.na(s) | s == 0, NA, value / s))
ridge <- c("vWP_lFDR", "RT_lFDR", "RE_lFDR", "RS_lasso", "RT_MPS")
mccA <- function(k) arr_of(k, "met")[, "MCC", ]

for (k in keys) {
  M <- mccA(k); I <- vec_of(k, "info")
  ## MB highest or near-highest
  best_other <- names(which.max(rowMeans(M[rownames(M) != "MB", ])))
  claim("S5.4-a", "Neighbourhood selection has the highest or near-highest MCC in every setting",
        k, paste("MB - best other (", best_other, ")"), M["MB", ] - M[best_other, ])
  ## vWP never worse than RT
  claim("S5.4-b", "The van Wieringen-Peeters version is never worse than the archetypal ridge",
        k, "vWP_lFDR - RT_lFDR", M["vWP_lFDR", ] - M["RT_lFDR", ])
  ## RE versus RT
  claim("S5.4-c", "The extrapolated estimator has a lower MCC than the plain ridge estimator in fourteen of the sixteen settings, by up to 0.07, and a higher one by less than 0.01 in the other two",
        k, "RE_lFDR - RT_lFDR", M["RE_lFDR", ] - M["RT_lFDR", ])
  ## ordering of the ridge layers
  o <- names(sort(rowMeans(M[ridge, ]), decreasing = TRUE))
  claim("S5.4-d", "Ranking of the ridge layers (best minus second best ridge layer)",
        k, paste(o[1], "-", o[2]), M[o[1], ] - M[o[2], ])
  claim("S5.4-e", "Second-best ridge layer minus third", k, paste(o[2], "-", o[3]),
        M[o[2], ] - M[o[3], ])
  ## lFDR false discovery proportions below 0.2
  A <- arr_of(k, "met")
  for (m in c("vWP_lFDR", "RT_lFDR", "RE_lFDR"))
    claim("S5.4-f", "Their [lFDR] false discovery proportions stay below 0.2", k,
          paste(m, "FDP"), A[m, "FDP", ])
  ## RT-MPS family-wise error
  fp <- round(A["RT_MPS", "FDP", ] * A["RT_MPS", "NE", ]) > 0
  claim("S5.4-g", "RT-MPS: share of replicates with at least one false edge", k,
        "FWER RT_MPS", fp, mean(fp), sqrt(mean(fp) * (1 - mean(fp)) / 100))
  ## rankings: RT best MCC against the better lasso path
  O <- vec_of(k, "orc")
  bl <- if (mean(O["GLASSO", ]) >= mean(O["MB", ])) "GLASSO" else "MB"
  claim("S5.5-a", "The ridge rankings attain an MCC within 0.01 of the best lasso-based path (band n >= 200, random n = 400)",
        k, paste("RT best MCC -", bl), O["RT", ] - O[bl, ])
  claim("S5.5-a2", "Same with the best ridge ranking (vWP)", k, paste("vWP best MCC -", bl),
        O["vWP", ] - O[bl, ])
  RK <- arr_of(k, "rk")
  claim("S5.5-b", "Exact separation share (RT)", k, "Sep.", RK["RT", "MARGIN", ] > 0,
        mean(RK["RT", "MARGIN", ] > 0),
        sqrt(mean(RK["RT", "MARGIN", ] > 0) * (1 - mean(RK["RT", "MARGIN", ] > 0)) / 100))
  claim("S5.5-c", "The areas under the ROC curve are above 0.95 for n >= 200", k,
        "RT AUC", RK["RT", "AUC", ])
  SE <- arr_of(k, "sens")
  claim("S5.5-d", "Best MCC of the ridge ranking, largest minus smallest grid penalty", k,
        "sens mcc(m=10) - mcc(m=0.3)", SE["mcc", 4, ] - SE["mcc", 1, ])
  claim("S5.5-e", "Grid penalty with the largest mean best MCC (band: peak between 1.2 and 4.5)",
        k, paste0("argmax at kappa = ", round(rowMeans(SE["kappa", , ])[which.max(rowMeans(SE["mcc", , ]))], 2)),
        NA, max(rowMeans(SE["mcc", , ])), NA)
  ## Section 5.3 leakage
  L <- arr_of(k, "leak")
  claim("S5.3-a", "For edges the observed and predicted means differ by at most 0.002", k,
        "obs - pred at distance 1", L["obs", "1", ] - L["pred", "1", ])
  claim("S5.3-b", "Population leakage at distance two between 0.014 and 0.055", k,
        "pred at distance 2", L["pred", "2", ])
  far <- if (all(is.na(L["obs", "3", ]))) L["obs", "4", ] else
    (L["obs", "3", ] + L["obs", "4", ]) / 2
  claim("S5.3-c", "Noise floor between 0.03 and 0.04 at distances of three or more", k,
        "obs at distance 3 and 4+ (mean of the two classes)", far)
  claim("S5.3-d", "Only for the band graph is the observed mean at distance two visibly above this floor",
        k, "obs(d=2) - obs(d=4+)", L["obs", "2", ] - L["obs", "4", ])
  for (dd in c("2", "4"))
    claim("S5.3-e", "The extrapolated estimator raises the observed means at almost all distances",
          k, paste0("rich - obs at distance ", dd), L["rich", dd, ] - L["obs", dd, ])
  ## Section 5.2
  claim("S5.2-a", "kappa_hat / lambda: 2.6-2.7 at n = 50, about 1.2 at n = 400", k,
        "kappa_hat / lambda_CV", I["kappa", ] / I["lambda", ])
  claim("S5.2-b", "Max deviation decreases from 0.18-0.19 (n = 50) to 0.13-0.14 (n = 400)",
        k, "max dev.", I["maxdev", ])
}
## Table 1 sentence: within 4% in every cell, below MC in all cells
for (i in seq_len(nrow(kg_pub))) {
  d <- kg[kg$lambda == kg_pub$lambda[i] & kg$n == kg_pub$n[i], ]
  dd <- (d$kappa_hat - d$kappa_mc_rep) / mean(d$kappa_mc_rep)
  claim("S5.2-c", "The estimate is within 4% of the Monte Carlo value in every cell",
        paste0("lambda=", kg_pub$lambda[i], " n=", kg_pub$n[i]),
        "(kappa_hat - kappa_mc) / kappa_mc, paired", dd)
}
claims <- do.call(rbind, cl_rows)
write.csv(claims, file.path(RES, "sim_claims.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Task 4: pair-level signed partial correlations, and the leakage summary
## ---------------------------------------------------------------------------
pairs <- do.call(rbind, lapply(keys, function(k) {
  g <- key_g(k); m <- models[[g]]
  ut <- upper.tri(m$E); ij <- which(ut, arr.ind = TRUE)
  H <- sapply(reps[[k]], `[[`, "rho_hat"); K <- sapply(reps[[k]], `[[`, "rho_kap")
  I <- vec_of(k, "info")
  data.frame(graph = g, n = key_n(k), i = ij[, 1], j = ij[, 2],
             distance = m$D[ut], edge = m$E[ut],
             rho_theta = pcor_of(m$Theta)[ut],
             rho_kappa_bar = pcor_of(solve(m$Sig + mean(I["kappa", ]) * diag(p_sim)))[ut],
             rho_kappa_mean = rowMeans(K),
             rho_hat_mean = rowMeans(H), rho_hat_se = apply(H, 1, se))
}))
write.csv(pairs, file.path(RES, "sim_pairs_signed.csv"), row.names = FALSE)
## how well the population ridge network and Theta predict the mean signed
## partial correlations: root mean square error by distance class, next to
## the mean Monte Carlo standard error of the pair means
pairs$dclass <- factor(pmin(pairs$distance, 4), levels = 1:4,
                       labels = c("1", "2", "3", "4+"))
pfit <- do.call(rbind, lapply(split(pairs, list(pairs$graph, pairs$n, pairs$dclass),
                                    drop = TRUE), function(d)
  data.frame(graph = d$graph[1], n = d$n[1], distance = d$dclass[1], pairs = nrow(d),
             rmse_vs_ridge_network = sqrt(mean((d$rho_hat_mean - d$rho_kappa_mean)^2)),
             rmse_vs_theta = sqrt(mean((d$rho_hat_mean - d$rho_theta)^2)),
             mean_mcse = mean(d$rho_hat_se),
             mean_rho_hat = mean(d$rho_hat_mean), mean_rho_kappa = mean(d$rho_kappa_mean),
             slope_on_ridge_network = if (sd(d$rho_kappa_mean) > 0)
               unname(coef(lm(rho_hat_mean ~ rho_kappa_mean, d))[2]) else NA)))
pfit <- pfit[order(match(pfit$graph, graphs), pfit$n, pfit$distance), ]
write.csv(pfit, file.path(RES, "sim_pairs_fit.csv"), row.names = FALSE)

leak <- do.call(rbind, lapply(keys, function(k) {
  L <- arr_of(k, "leak")
  data.frame(graph = key_g(k), n = key_n(k), distance = c("1", "2", "3", "4+"),
             observed = apply(L["obs", , ], 1, mean), se_observed = apply(L["obs", , ], 1, se),
             predicted = apply(L["pred", , ], 1, mean),
             extrapolated = apply(L["rich", , ], 1, mean))
}))
write.csv(leak, file.path(RES, "sim_leakage.csv"), row.names = FALSE)

## extrapolated matrix: how often it is not positive definite
pdr <- do.call(rbind, lapply(keys, function(k) {
  P <- vec_of(k, "pd_rich")
  data.frame(graph = key_g(k), n = key_n(k), share_not_pd = mean(P["min_eig", ] <= 0),
             share_nonpos_diag = mean(P["min_diag", ] <= 0))
}))
write.csv(pdr, file.path(RES, "sim_extrap_pd.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Task 5: p >> n block
## ---------------------------------------------------------------------------
hd_summary <- function(pp) {
  fs <- list.files(WORK, pattern = sprintf("^hd_p%d_.*\\.rds$", pp), full.names = TRUE)
  if (!length(fs)) return(NULL)
  info <- do.call(rbind, lapply(fs, function(f) {
    b <- sub("\\.rds$", "", basename(f)); parts <- strsplit(b, "_")[[1]]
    data.frame(file = f, graph = parts[3], n = as.integer(parts[4]))
  }))
  do.call(rbind, lapply(split(info, list(info$graph, info$n), drop = TRUE), function(d) {
    R <- lapply(d$file, readRDS)
    ok <- !sapply(R, function(x) !is.null(x$error))
    R <- R[ok]
    cols <- function(q) do.call(cbind, lapply(R, `[[`, q))
    I <- cols("info"); O <- cols("orc"); RK <- cols("rk"); P <- cols("pop")
    TM <- cols("seconds"); rownames(TM) <- sub("[.].*$", "", rownames(TM))
    A <- simplify2array(lapply(R, `[[`, "met"), higher = TRUE)
    if (length(dim(A)) == 2) A <- array(A, c(dim(A), 1), c(dimnames(A), list(NULL)))
    fp <- matrix(round(A[, "FDP", ] * A[, "NE", ]), nrow = dim(A)[1],
                 dimnames = list(dimnames(A)[[1]], NULL))
    base <- data.frame(p = pp, graph = d$graph[1], n = d$n[1], reps_ok = sum(ok),
      reps_failed = sum(!ok),
      edges = mean(I["n_edges", ]), pc_min = mean(I["pc_min", ]),
      pc_median = mean(I["pc_med", ]), pc_max = mean(I["pc_max", ]),
      cond = median(I["cond", ]),
      lambda_cv = mean(I["lambda", ]), kappa_hat = mean(I["kappa", ]),
      sd_kappa_hat = sd(I["kappa", ]), kappa_mc = mean(I["kappa_mc", ]),
      rel_diff_kappa = mean((I["kappa", ] - I["kappa_mc", ]) / I["kappa_mc", ]),
      se_rel_diff = se((I["kappa", ] - I["kappa_mc", ]) / I["kappa_mc", ]),
      floor = mean(I["floor", ]), floor_kyfan = mean(I["floor_kyfan", ]),
      zeta_pop = mean(P["zeta", ]), ell_pop = mean(P["ell", ]),
      share_ell_gt_zeta = mean(P["ell", ] > P["zeta", ]),
      best_MCC_pop = mean(P["mcc", ]),
      best_MCC_RT = mean(O["RT", ]), se_best_MCC_RT = se(O["RT", ]),
      AUC_RT = mean(RK["AUC", ]), sep_RT = mean(RK["MARGIN", ] > 0),
      best_MCC_GLASSO = mean(O["GLASSO", ]), best_MCC_MB = mean(O["MB", ]),
      sec_ridge = mean(TM["ridge", ]), sec_mps = mean(TM["mps", ]),
      sec_glasso = mean(TM["glasso", ]), sec_mb = mean(TM["mb", ]))
    for (m in dimnames(A)[[1]]) {
      base[[paste0(m, "_TPR")]] <- mean(A[m, "TPR", ])
      base[[paste0(m, "_FDP")]] <- mean(A[m, "FDP", ])
      base[[paste0(m, "_MCC")]] <- mean(A[m, "MCC", ])
      base[[paste0(m, "_seMCC")]] <- se(A[m, "MCC", ])
      base[[paste0(m, "_NE")]] <- mean(A[m, "NE", ])
      base[[paste0(m, "_FWER")]] <- mean(fp[m, ] > 0)
    }
    base
  }))
}
## row-bind tables whose columns differ (p = 1000 has the ridge parts only)
rbind_fill <- function(a, b) {
  if (is.null(a)) return(b); if (is.null(b)) return(a)
  for (v in setdiff(names(a), names(b))) b[[v]] <- NA
  for (v in setdiff(names(b), names(a))) a[[v]] <- NA
  rbind(a, b[names(a)])
}
hd <- rbind_fill(hd_summary(400), hd_summary(1000))
if (!is.null(hd)) {
  hd <- hd[order(hd$p, match(hd$graph, c("band", "hub", "random")), hd$n), ]
  write.csv(hd, file.path(RES, "sim_highdim.csv"), row.names = FALSE)
  ## add the p >> n error rates to sim_fwer.csv
  for (m in c("RT_lFDR", "RT_MPS", "GLASSO", "MB")) {
    if (is.null(hd[[paste0(m, "_TPR")]])) next
    fwer <- rbind(fwer, data.frame(design = "high-dimensional", p = hd$p,
      graph = hd$graph, n = hd$n, method = m, TPR = hd[[paste0(m, "_TPR")]],
      FDP = hd[[paste0(m, "_FDP")]], MCC = hd[[paste0(m, "_MCC")]],
      seMCC = hd[[paste0(m, "_seMCC")]], NE = hd[[paste0(m, "_NE")]],
      FWER = hd[[paste0(m, "_FWER")]],
      seFWER = sqrt(hd[[paste0(m, "_FWER")]] * (1 - hd[[paste0(m, "_FWER")]]) / hd$reps_ok),
      mean_FP = NA))
  }
  fwer <- fwer[!(fwer$design == "high-dimensional" & is.na(fwer$TPR)), ]
}
write.csv(fwer, file.path(RES, "sim_fwer.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo_summaries.txt"))
cat("summaries written\n")
