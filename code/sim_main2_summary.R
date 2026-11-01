###############################################################################
## Tables for the second extension of the main simulation design
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Purpose: reads the replicates written by sim_extensions.R
##          (../results/ext_work/, replicates 1-100) and sim_main2.R
##          (../results/ext_work3/), and writes ../results/sim_main2_*.csv.
##          Means are given with Monte Carlo standard errors; tables are made
##          for the first 100 replicates (the published ones) and for all 500.
## Usage  : Rscript sim_main2_summary.R      (from code/)
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages({ library(huge); library(igraph) })
OUT  <- ".."
RES  <- file.path(OUT, "results")
W0   <- file.path(RES, "ext_work")
W3   <- file.path(RES, "ext_work3")
for (e in parse("sim_ridge_support.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function")))
    eval(e)
}
p_sim  <- 100
n_grid <- c(50, 100, 200, 400)
graphs <- c("band", "random", "hub", "scale-free")
NREP <- 500; NOLD <- 100
se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
## best MCC over all thresholds of |s| (s on the upper triangle, e edges)
oracle_mcc_threshold_v <- function(s, e) {
  o <- order(abs(s), decreasing = TRUE)
  tp <- cumsum(e[o]); fp <- cumsum(!e[o])
  P <- sum(e); N <- sum(!e); fn <- P - tp; tn <- N - fp
  den <- sqrt(tp + fp) * sqrt(tp + fn) * sqrt(tn + fp) * sqrt(tn + fn)
  max(ifelse(den > 0, (tp * tn - fp * fn) / den, 0))
}

## ---------------------------------------------------------------------------
## 1. Read the replicates
## ---------------------------------------------------------------------------
keys <- as.vector(outer(graphs, n_grid, paste, sep = "_"))
keys <- keys[order(match(sub("_.*", "", keys), graphs),
                   as.integer(sub(".*_", "", keys)))]
key_g <- function(k) sub("_[0-9]+$", "", k)
key_n <- function(k) as.integer(sub(".*_", "", k))
A <- list(); X0 <- list()
for (k in keys) {
  A[[k]] <- list(); X0[[k]] <- list()
  for (r in seq_len(NREP)) {
    f <- file.path(W3, sprintf("main2_%s_rep%03d.rds", k, r))
    if (!file.exists(f)) next
    a <- readRDS(f)
    if (r <= NOLD) a$ext <- readRDS(file.path(W0, sprintf("main_%s_rep%03d.rds",
                                                         k, r)))$res
    A[[k]][[r]] <- a
  }
  A[[k]] <- A[[k]][!vapply(A[[k]], is.null, logical(1))]
}
nrep_of <- sapply(A, length)
cat("replicates per setting:\n"); print(nrep_of)

## ---------------------------------------------------------------------------
## 2. Checks: rebuilt replicates against the stored ones
## ---------------------------------------------------------------------------
chk <- do.call(rbind, lapply(keys, function(k) {
  C <- sapply(A[[k]], `[[`, "chk")
  data.frame(graph = key_g(k), n = key_n(k), reps = ncol(C),
             t(apply(abs(C), 1, max)), row.names = NULL)
}))
write.csv(chk, file.path(RES, "sim_main2_check.csv"), row.names = FALSE)

## helper: one statistic per replicate as a matrix (rows = replicates)
grab <- function(k, f, reps = NREP) {
  L <- A[[k]][seq_len(min(reps, length(A[[k]])))]
  L <- lapply(L, f)
  nm <- unique(unlist(lapply(L, names)))
  t(sapply(L, function(v) { out <- setNames(rep(NA_real_, length(nm)), nm)
                            out[names(v)] <- v; out }))
}
summ_block <- function(M, k, reps, prefix = "") {
  data.frame(graph = key_g(k), n = key_n(k), reps = reps,
             method = colnames(M), mean = colMeans(M, na.rm = TRUE),
             se = apply(M, 2, se), nobs = colSums(!is.na(M)), row.names = NULL)
}

## ---------------------------------------------------------------------------
## 3. Selection: TPR, FDP, MCC, size, FWER
## ---------------------------------------------------------------------------
sel_all <- function(a) rbind(a$ext$met, a$ext$met2, a$met3)
selection <- do.call(rbind, lapply(c(NOLD, NREP), function(R) do.call(rbind,
  lapply(keys, function(k) {
    L <- A[[k]][seq_len(min(R, length(A[[k]])))]
    ms <- unique(unlist(lapply(L, function(a) rownames(sel_all(a)))))
    ms <- setdiff(ms, "vWP_lFDR_chk")
    do.call(rbind, lapply(ms, function(m) {
      M <- t(sapply(L, function(a) { s <- sel_all(a)
        if (m %in% rownames(s)) s[m, ] else rep(NA, 4) }))
      ok <- !is.na(M[, 1])
      M <- M[ok, , drop = FALSE]
      data.frame(graph = key_g(k), n = key_n(k), reps = R, method = m,
                 nobs = nrow(M), TPR = mean(M[, 1]), FDP = mean(M[, 2]),
                 seFDP = se(M[, 2]), MCC = mean(M[, 3]), seMCC = se(M[, 3]),
                 NE = mean(M[, 4]), seNE = se(M[, 4]),
                 FWER = mean(M[, 2] > 0), row.names = NULL)
    }))
  }))))
write.csv(selection, file.path(RES, "sim_main2_selection.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4. Rankings: best MCC, matched size, AUC
## ---------------------------------------------------------------------------
ranking <- do.call(rbind, lapply(c(NOLD, NREP), function(R) do.call(rbind,
  lapply(keys, function(k) {
    L <- A[[k]][seq_len(min(R, length(A[[k]])))]
    rk <- lapply(L, `[[`, "rank"); pa <- lapply(L, `[[`, "paths")
    ms <- unique(unlist(lapply(rk, rownames)))
    one <- function(m, get, cols) {
      M <- do.call(rbind, lapply(seq_along(L), function(i) {
        x <- get(i, m); if (is.null(x)) rep(NA, length(cols)) else x[cols] }))
      M <- M[!is.na(M[, 1]), , drop = FALSE]
      data.frame(graph = key_g(k), n = key_n(k), reps = R, method = m,
                 nobs = nrow(M), t(colMeans(M, na.rm = TRUE)),
                 t(setNames(apply(M, 2, se), paste0("se_", cols))),
                 row.names = NULL)
    }
    cols <- c("best", "mcc_m", "tpr_m", "mcc_liu", "tpr_liu", "fdp_liu")
    r1 <- do.call(rbind, lapply(ms, one, get = function(i, m)
      if (m %in% rownames(rk[[i]])) rk[[i]][m, ] else NULL, cols = c(cols, "auc")))
    pcols <- c("best", "mcc_m.MCC", "tpr_m.TPR", "mcc_liu.MCC", "tpr_liu.TPR",
               "fdp_liu.FDP")
    r2 <- do.call(rbind, lapply(c("GLASSO", "MB"), function(m) {
      d <- one(m, function(i, m) pa[[i]][m, ], c(pcols, "size_m", "size_liu"))
      names(d) <- sub("\\.(MCC|TPR|FDP)", "", names(d)); d }))
    r1$size_m <- NA; r1$size_liu <- NA; r1$se_size_m <- NA; r1$se_size_liu <- NA
    r2$auc <- NA; r2$se_auc <- NA
    rbind(r1, r2[, names(r1)])
  }))))
write.csv(ranking, file.path(RES, "sim_main2_ranking.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. FDR rules at the nominal level 0.1
## ---------------------------------------------------------------------------
fdr_methods <- c("HS_BH10", "Liu2013", "RT_Fdr10", "RT_BH10", "vWP_Fdr10",
                 "vWPt_Fdr10", "GeneNet_BH10", "GeneNet_Fdr10", "beam_BH10")
fdr <- subset(selection, method %in% fdr_methods)
fdr$excess_z <- (fdr$FDP - 0.1) / fdr$seFDP
write.csv(fdr, file.path(RES, "sim_main2_fdr.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 6. Tuning: penalties and effective penalties of all ridge-type methods
## ---------------------------------------------------------------------------
tuning <- do.call(rbind, lapply(keys, function(k) {
  Tm <- grab(k, function(a) a$tune)
  d <- summ_block(Tm, k, nrow(Tm))
  ## shares that describe the choices
  extra <- data.frame(graph = key_g(k), n = key_n(k), reps = nrow(Tm),
    method = c("share_cv100_differs", "share_omega_stab_30", "share_omega_stabp_30",
               "share_vWPt_at_edge"),
    mean = c(mean(abs(Tm[, "lambda_cv100"] / Tm[, "lambda"] - 1) > 0.02),
             mean(Tm[, "omega_stab"] == 30), mean(Tm[, "omega_stabp"] == 30),
             mean(Tm[, "vWPt_at_edge"])),
    se = NA, nobs = nrow(Tm))
  rbind(d, extra)
}))
write.csv(tuning, file.path(RES, "sim_main2_tuning.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 7. Paired differences with multiplicity control over the 16 settings
## ---------------------------------------------------------------------------
## per-replicate value of a selection MCC or a ranking summary
v_sel <- function(m, col = "MCC") function(a) { s <- sel_all(a)
  if (m %in% rownames(s)) s[m, col] else NA }
v_rank <- function(m, col = "best") function(a) {
  if (m %in% rownames(a$rank)) a$rank[m, col]
  else if (m %in% rownames(a$paths)) a$paths[m, if (col == "best") "best" else
    paste0(col, ".", toupper(sub("_.*", "", col)))] else NA }
ridge_rules <- c("vWP_lFDR", "RT_lFDR", "RE_lFDR", "RS_lasso", "RT_MPS")
sel_main <- c(ridge_rules, "HS_BH10", "Liu2013", "GLASSO", "MB")
## comparisons: name, value of A, value of B (B may be "best of a set")
cmp <- list(
  list("MB - best other rule (MCC)", v_sel("MB"), setdiff(sel_main, "MB")),
  list("vWP-lFDR - RT-lFDR (MCC)", v_sel("vWP_lFDR"), "RT_lFDR"),
  list("RE-lFDR - RT-lFDR (MCC)", v_sel("RE_lFDR"), "RT_lFDR"),
  list("vWP-lFDR 0.9 - 0.99 (MCC)", v_sel("vWP_lFDR"), "vWP_lFDR99"),
  list("HS - best ridge rule (MCC)", v_sel("HS_BH10"), ridge_rules),
  list("Liu - best ridge rule (MCC)", v_sel("Liu2013"), ridge_rules),
  list("Liu - MB (MCC)", v_sel("Liu2013"), "MB"),
  list("RS - MB (MCC)", v_sel("RS_lasso"), "MB"),
  list("vWP-target - vWP-null, lFDR 0.9 (MCC)", v_sel("vWPt_lFDR90"), "vWP_lFDR"),
  list("GeneNet lFDR 0.9 - vWP-lFDR (MCC)", v_sel("GeneNet_lFDR90"), "vWP_lFDR"),
  list("GeneNet lFDR 0.8 - vWP-lFDR (MCC)", v_sel("GeneNet_lFDR80"), "vWP_lFDR"),
  list("beam BH 0.1 - MB (MCC)", v_sel("beam_BH10"), "MB"),
  list("beam BH 0.1 - Liu (MCC)", v_sel("beam_BH10"), "Liu2013"),
  list("Ranking: Marginal - RT (best MCC)", v_rank("Marginal"), "rank:RT"),
  list("Ranking: MB path - RT (best MCC)", v_rank("MB"), "rank:RT"),
  list("Ranking: GLASSO path - RT (best MCC)", v_rank("GLASSO"), "rank:RT"),
  list("Ranking: vWP - RT (best MCC)", v_rank("vWP"), "rank:RT"),
  list("Ranking: RE - RT (best MCC)", v_rank("RE"), "rank:RT"),
  list("Ranking: vWP-target - vWP (best MCC)", v_rank("vWPt"), "rank:vWP"),
  list("Ranking: GeneNet - RT (best MCC)", v_rank("GeneNet"), "rank:RT"),
  list("Ranking: beam - RT (best MCC)", v_rank("beam"), "rank:RT"),
  list("Ranking: MB path - GeneNet (best MCC)", v_rank("MB"), "rank:GeneNet"),
  list("Ranking: MB path - beam (best MCC)", v_rank("MB"), "rank:beam"),
  list("Ranking: MB - RT at true size (MCC)", v_rank("MB", "mcc_m"), "rankm:RT"),
  list("Ranking: oracle leakage-only correction - RT (best MCC)", v_rank("ORC_leak"), "rank:RT"),
  list("Ranking: oracle full correction - RT (best MCC)", v_rank("ORC"), "rank:RT"),
  list("Ranking: RT at CV 100-grid - RT (best MCC)", v_rank("RT_cv100"), "rank:RT"),
  list("Ranking: RT at stability penalty - RT (best MCC)", v_rank("RT_stab"), "rank:RT"),
  list("Ranking: RT at 10 x CV - RT (best MCC)", v_rank("RT_x10"), "rank:RT"),
  list("Ranking: MB path - RT at stability penalty (best MCC)", v_rank("MB"), "rank:RT_stab"))
b_value <- function(spec, a, k_means) {
  if (length(spec) > 1) return(v_sel(spec[which.max(k_means[spec])])(a))
  if (grepl("^rankm:", spec)) return(v_rank(sub("^rankm:", "", spec), "mcc_m")(a))
  if (grepl("^rank:", spec)) return(v_rank(sub("^rank:", "", spec))(a))
  v_sel(spec)(a)
}
zc <- qnorm(1 - 0.05 / (2 * length(keys)))     # simultaneous over 16 settings
mult <- do.call(rbind, lapply(c(NOLD, NREP), function(R) do.call(rbind,
  lapply(cmp, function(cm) {
    d <- do.call(rbind, lapply(keys, function(k) {
      L <- A[[k]][seq_len(min(R, length(A[[k]])))]
      ## means of the rules, for "best of a set" (chosen per setting)
      km <- sapply(unique(c(sel_main, cm[[3]][!grepl(":", cm[[3]])])), function(m)
        mean(sapply(L, v_sel(m)), na.rm = TRUE))
      bname <- if (length(cm[[3]]) > 1) names(which.max(km[cm[[3]]])) else cm[[3]]
      x <- sapply(L, cm[[2]]) - sapply(L, function(a) b_value(cm[[3]], a, km))
      x <- x[!is.na(x)]
      m <- mean(x); s <- se(x)
      data.frame(comparison = cm[[1]], reps = R, graph = key_g(k), n = key_n(k),
                 B = bname, nobs = length(x), diff = m, se = s, z = m / s,
                 p = 2 * pnorm(-abs(m / s)), lo_sim = m - zc * s, hi_sim = m + zc * s)
    }))
    d$p[!is.finite(d$p)] <- 1
    d$p_holm <- p.adjust(d$p, "holm")
    d$verdict <- ifelse(abs(d$diff) <= 2 * d$se, "tie (within 2 MCSE)",
                 ifelse(d$p_holm < 0.05, ifelse(d$diff > 0, "higher (Holm)", "lower (Holm)"),
                        ifelse(d$diff > 0, "higher (unadjusted only)",
                               "lower (unadjusted only)")))
    d
  }))))
write.csv(mult, file.path(RES, "sim_main2_multiplicity.csv"), row.names = FALSE)
mult_count <- aggregate(cbind(higher = verdict == "higher (Holm)",
                              lower = verdict == "lower (Holm)",
                              tie = verdict == "tie (within 2 MCSE)",
                              unadj = grepl("unadjusted", verdict),
                              raw_sig = p < 0.05) ~ comparison + reps,
                        data = mult, FUN = sum)
mult_count <- mult_count[order(match(mult_count$comparison,
                                     sapply(cmp, `[[`, 1)), mult_count$reps), ]
write.csv(mult_count, file.path(RES, "sim_main2_multiplicity_counts.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## 8. Leakage at the sample level and the extrapolation diagnostic
## ---------------------------------------------------------------------------
pair_diff <- function(x, y) { d <- x - y; c(diff = mean(d), se = se(d)) }
leak_tab <- do.call(rbind, lapply(keys, function(k) {
  L <- A[[k]]
  b <- function(m, col = "best") sapply(L, function(a) a$rank[m, col])
  pr <- t(sapply(L, function(a) a$pop_re["om0", ]))
  data.frame(graph = key_g(k), n = key_n(k), reps = length(L),
             RT = mean(b("RT")), ORC_leak = mean(b("ORC_leak")),
             ORC = mean(b("ORC")),
             d_leak = pair_diff(b("ORC_leak"), b("RT"))[1],
             se_d_leak = pair_diff(b("ORC_leak"), b("RT"))[2],
             d_full = pair_diff(b("ORC"), b("RT"))[1],
             se_d_full = pair_diff(b("ORC"), b("RT"))[2],
             RT_m = mean(b("RT", "mcc_m")), ORC_leak_m = mean(b("ORC_leak", "mcc_m")),
             ORC_m = mean(b("ORC", "mcc_m")),
             pop_margin = mean(pr[, "zeta_RT"] - pr[, "ell_RT"]),
             row.names = NULL)
}))
## the p = 400 block
hd_graphs <- c("band", "hub", "random")
HD <- list()
for (g in hd_graphs) for (n in c(50, 100)) {
  k <- paste(g, n, sep = "_")
  fs <- file.path(W3, sprintf("hd2_p400_%s_%d_rep%03d.rds", g, n, 1:100))
  HD[[k]] <- lapply(fs[file.exists(fs)], readRDS)
}
leak_hd <- do.call(rbind, lapply(names(HD), function(k) {
  L <- HD[[k]]
  b <- function(m, col = "best") sapply(L, function(a) a$rank[m, col])
  data.frame(graph = paste0(key_g(k), ", p = 400"), n = key_n(k), reps = length(L),
             RT = mean(b("RT")), ORC_leak = mean(b("ORC_leak")), ORC = mean(b("ORC")),
             d_leak = pair_diff(b("ORC_leak"), b("RT"))[1],
             se_d_leak = pair_diff(b("ORC_leak"), b("RT"))[2],
             d_full = pair_diff(b("ORC"), b("RT"))[1],
             se_d_full = pair_diff(b("ORC"), b("RT"))[2],
             RT_m = mean(b("RT", "mcc_m")), ORC_leak_m = mean(b("ORC_leak", "mcc_m")),
             ORC_m = mean(b("ORC", "mcc_m")),
             pop_margin = mean(sapply(L, function(a) a$pop["RT", "zeta"] - a$pop["RT", "ell"])),
             row.names = NULL)
}))
write.csv(rbind(leak_tab, leak_hd), file.path(RES, "sim_main2_leakage_oracle.csv"),
          row.names = FALSE)

## extrapolation: main design at lambda_CV and at the larger penalties
ex_main <- do.call(rbind, lapply(keys, function(k) do.call(rbind,
  lapply(c("om0", "om1", "om3", "om10"), function(om) {
    P <- t(sapply(A[[k]], function(a) a$pop_re[om, ]))
    data.frame(graph = key_g(k), n = key_n(k), reps = nrow(P),
               lambda1 = if (om == "om0") "CV" else sub("om", "", om),
               kappa1 = mean(P[, "kappa1"]), ratio_l2_l1 = mean(P[, "ratio"]),
               max_ratio = max(P[, "ratio"]),
               RT = mean(P[, "RT_best"]), RE = mean(P[, "RE_best"]),
               d_RE_RT = mean(P[, "RE_best"] - P[, "RT_best"]),
               se_d = se(P[, "RE_best"] - P[, "RT_best"]),
               pop_margin_RT = mean(P[, "zeta_RT"] - P[, "ell_RT"]),
               pop_margin_RE = mean(P[, "zeta_RE"] - P[, "ell_RE"]),
               share_neg_margin_RT = mean(P[, "zeta_RT"] < P[, "ell_RT"]),
               share_neg_margin_RE = mean(P[, "zeta_RE"] < P[, "ell_RE"]),
               pop_mcc_RT = mean(P[, "mcc_RT"]), pop_mcc_RE = mean(P[, "mcc_RE"]),
               row.names = NULL)
  }))))
write.csv(ex_main, file.path(RES, "sim_main2_extrap_main.csv"), row.names = FALSE)

## root-searched against grid second penalty, at lambda_CV
ex_root <- rbind(
  do.call(rbind, lapply(keys, function(k) {
    I <- t(sapply(A[[k]], `[[`, "re_info"))
    b <- function(m) sapply(A[[k]], function(a) a$rank[m, "best"])
    data.frame(graph = key_g(k), n = key_n(k), reps = nrow(I),
               ratio_mean = mean(I[, "ratio"]), ratio_min = min(I[, "ratio"]),
               ratio_max = max(I[, "ratio"]),
               share_pd_guarantee = mean(I[, "pd_guarantee"]),
               share_not_pd = mean(I[, "min_eig"] <= 0),
               h_g30_min = min(I[, "h_g30"]), h_g30_max = max(I[, "h_g30"]),
               RT = mean(b("RT")), RE_root = mean(b("RE")), RE_g30 = mean(b("RE_g30")),
               RE_g100 = mean(b("RE_g100")),
               d_root_RT = mean(b("RE") - b("RT")), se_root_RT = se(b("RE") - b("RT")),
               d_g30_root = mean(b("RE_g30") - b("RE")),
               se_g30_root = se(b("RE_g30") - b("RE")),
               maxabs_g30_root = max(abs(b("RE_g30") - b("RE"))), row.names = NULL)
  })),
  do.call(rbind, lapply(names(HD), function(k) {
    L <- HD[[k]]
    I <- t(sapply(L, `[[`, "re_info"))
    b <- function(m) sapply(L, function(a) a$rank[m, "best"])
    data.frame(graph = paste0(key_g(k), ", p = 400"), n = key_n(k), reps = nrow(I),
               ratio_mean = mean(I[, "ratio"]), ratio_min = min(I[, "ratio"]),
               ratio_max = max(I[, "ratio"]),
               share_pd_guarantee = mean(I[, "pd_guarantee"]),
               share_not_pd = mean(I[, "min_eig"] <= 0),
               h_g30_min = min(I[, "h_g30"]), h_g30_max = max(I[, "h_g30"]),
               RT = mean(b("RT")), RE_root = mean(b("RE")), RE_g30 = mean(b("RE_g30")),
               RE_g100 = NA,
               d_root_RT = mean(b("RE") - b("RT")), se_root_RT = se(b("RE") - b("RT")),
               d_g30_root = mean(b("RE_g30") - b("RE")),
               se_g30_root = se(b("RE_g30") - b("RE")),
               maxabs_g30_root = max(abs(b("RE_g30") - b("RE"))), row.names = NULL)
  })))
write.csv(ex_root, file.path(RES, "sim_main2_extrap_root.csv"), row.names = FALSE)

## population margin of the extrapolated target in the p = 400 block
ex_hd_pop <- do.call(rbind, lapply(names(HD), function(k) {
  L <- HD[[k]]
  P <- t(sapply(L, function(a) c(a$pop["RT", ], a$pop["RE", ])))
  colnames(P) <- c("zeta_RT", "ell_RT", "mcc_RT", "zeta_RE", "ell_RE", "mcc_RE")
  data.frame(graph = paste0(key_g(k), ", p = 400"), n = key_n(k), reps = nrow(P),
             pop_margin_RT = mean(P[, 1] - P[, 2]), pop_margin_RE = mean(P[, 4] - P[, 5]),
             share_neg_RT = mean(P[, 1] < P[, 2]), share_neg_RE = mean(P[, 4] < P[, 5]),
             pop_mcc_RT = mean(P[, 3]), pop_mcc_RE = mean(P[, 6]), row.names = NULL)
}))
write.csv(ex_hd_pop, file.path(RES, "sim_main2_extrap_hd_pop.csv"), row.names = FALSE)

## cross-validation in the p = 400 block: choices and the curve
hd_cv <- do.call(rbind, lapply(names(HD), function(k) {
  L <- HD[[k]]
  lam30 <- sapply(L, function(a) a$re_info[["lambda1"]])
  lam100 <- sapply(L, `[[`, "lambda_cv100")
  C30 <- sapply(L, `[[`, "cv30")
  grid30 <- exp(seq(log(1e-2), log(20), length.out = 30))
  data.frame(graph = key_g(k), n = key_n(k), reps = length(L),
             lambda_cv30_values = paste(sort(unique(signif(lam30, 4))), collapse = " "),
             lambda_cv100_mean = mean(lam100), lambda_cv100_min = min(lam100),
             lambda_cv100_max = max(lam100),
             kappa_cv100 = mean(sapply(L, `[[`, "kappa_cv100")),
             RT_cv100_best = mean(sapply(L, function(a) a$rank["RT_cv100", "best"])),
             RT_best = mean(sapply(L, function(a) a$rank["RT", "best"])),
             curve_rel_at_neighbours = paste(signif(
               rowMeans(sweep(C30, 2, apply(C30, 2, min), "/"))[16:20], 4), collapse = " "),
             row.names = NULL)
}))
write.csv(hd_cv, file.path(RES, "sim_main2_hd_cv.csv"), row.names = FALSE)
hd_curve <- do.call(rbind, lapply(names(HD), function(k) {
  C <- sapply(HD[[k]], `[[`, "cv100")
  data.frame(graph = key_g(k), n = key_n(k),
             lambda = exp(seq(log(1e-2), log(20), length.out = 100)),
             mean_cv_loss = rowMeans(C), mean_rel_to_min = rowMeans(sweep(C, 2,
               apply(C, 2, min), "/")))
}))
write.csv(hd_curve, file.path(RES, "sim_main2_hd_cv_curve.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 9. Population ridge network for other ridge-type estimators
## ---------------------------------------------------------------------------
## For each estimator the replicate mean of the signed partial correlations
## is compared with the partial correlations of (Sigma + kappa' I)^{-1}; kappa'
## minimizes the root mean square difference over all pairs.
## the random and scale-free generators use random numbers: rebuild them with
## the generator of the simulation, and check the number of edges
RNGkind("L'Ecuyer-CMRG")
models <- lapply(setNames(graphs, graphs), make_model, p = p_sim)
for (k in keys) {
  m_true <- A[[k]][[1]]$tune[["m_true"]]
  E <- models[[key_g(k)]]$E
  stopifnot(sum(E[upper.tri(E)]) == m_true)
}
bridge <- do.call(rbind, lapply(keys, function(k) {
  m <- models[[key_g(k)]]; ut <- upper.tri(m$E)
  dcl <- pmin(m$D[ut], 4)
  L <- A[[k]][seq_len(min(NOLD, length(A[[k]])))]
  est <- unique(unlist(lapply(L, function(a) names(a$signed))))
  Tm <- grab(k, function(a) a$tune, NOLD)
  do.call(rbind, lapply(est, function(e) {
    H <- sapply(L, function(a) a$signed[[e]] / 1e5)
    if (is.list(H)) H <- do.call(cbind, H[!vapply(H, is.null, logical(1))])
    h <- rowMeans(H); hse <- apply(H, 1, sd) / sqrt(ncol(H))
    pop <- function(kp) pcor_of(solve(m$Sig + kp * diag(p_sim)))[ut]
    rms <- function(kp) sqrt(mean((h - pop(kp))^2))
    o <- optimize(function(lk) rms(exp(lk)), c(log(1e-3), log(1e3)))
    kp <- exp(o$minimum); d <- h - pop(kp)
    ## a distance class may be empty (the hub graph has no pairs at distance 3)
    by_d <- tapply(d, factor(dcl, levels = 1:4), function(x) sqrt(mean(x^2)))
    se_d <- tapply(hse, factor(dcl, levels = 1:4), mean)
    ref <- switch(e, RT = mean(Tm[, "kappa"]), GeneNet = mean(Tm[, "kappa_GeneNet"]),
                  beam = mean(Tm[, "kappa_beam"]), NA)
    rms_ref <- if (is.na(ref)) NA else rms(ref)
    ## the paper's comparison for RT: mean of rho^(kappa_hat) over replicates
    rms_paper <- if (e == "RT") sqrt(mean((h - rowMeans(sapply(L, function(a)
      a$ext$rho_kap)))^2)) else NA
    data.frame(graph = key_g(k), n = key_n(k), estimator = e, reps = ncol(H),
               kappa_fit = kp, rms_fit = o$objective, mean_mcse = mean(hse),
               rms_d1 = by_d[["1"]], rms_d2 = by_d[["2"]], rms_d3 = by_d[["3"]],
               rms_d4 = by_d[["4"]], mcse_d1 = se_d[["1"]], mcse_d2 = se_d[["2"]],
               kappa_hat_mean = ref, rms_at_kappa_hat = rms_ref,
               rms_paper_RT = rms_paper,
               rms_vs_theta = sqrt(mean((h - pcor_of(m$Theta)[ut])^2)),
               row.names = NULL)
  }))
}))
write.csv(bridge, file.path(RES, "sim_main2_bridge.csv"), row.names = FALSE)

## Oracle ranking on a common scale (replicates 1-100). The correction
## rho_hat - {rho^(kappa_hat) - P(Theta)} puts the edges back on the scale of
## Theta while the noise stays on the shrunken ridge scale. Here P(Theta) is
## multiplied by the mean shrinkage of the edges, sbar = mean over edges of
## |rho^(kappa_hat)| / |P(Theta)|, so that signal and noise keep the ridge
## scale: the ranking removes the leakage and the unequal shrinkage of the
## edges, and nothing else.
orc_scaled <- do.call(rbind, lapply(keys, function(k) {
  m <- models[[key_g(k)]]; ut <- upper.tri(m$E); e <- m$E[ut]
  Pt <- pcor_of(m$Theta)[ut]
  L <- A[[k]][seq_len(min(NOLD, length(A[[k]])))]
  v <- t(sapply(L, function(a) {
    rh <- a$signed$RT / 1e5
    rk <- pcor_of(solve(m$Sig + a$tune[["kappa"]] * diag(p_sim)))[ut]
    sbar <- mean(abs(rk[e]) / abs(Pt[e]))
    c(RT = oracle_mcc_threshold_v(rh, e),
      ORC_scaled = oracle_mcc_threshold_v(rh - rk + sbar * Pt, e),
      ORC_leak = a$rank["ORC_leak", "best"], ORC = a$rank["ORC", "best"],
      sbar = sbar)
  }))
  data.frame(graph = key_g(k), n = key_n(k), reps = nrow(v), t(colMeans(v)),
             d_scaled = mean(v[, "ORC_scaled"] - v[, "RT"]),
             se_d_scaled = se(v[, "ORC_scaled"] - v[, "RT"]), row.names = NULL)
}))
write.csv(orc_scaled, file.path(RES, "sim_main2_leakage_oracle_scaled.csv"),
          row.names = FALSE)

## ---------------------------------------------------------------------------
## 10. Table 2 of the paper: relative error with more digits
## ---------------------------------------------------------------------------
kg <- read.csv(file.path(RES, "sim_kappa_grid_reps.csv"))
kgt <- do.call(rbind, lapply(split(kg, list(kg$lambda, kg$n)), function(d)
  data.frame(lambda = d$lambda[1], n = d$n[1], kappa_mc = mean(d$kappa_mc_rep),
             kappa_hat = mean(d$kappa_hat),
             rel_err_pct = 100 * (mean(d$kappa_hat) / mean(d$kappa_mc_rep) - 1),
             diff = mean(d$kappa_hat - d$kappa_mc_rep),
             se_diff = se(d$kappa_hat - d$kappa_mc_rep))))
kgt <- kgt[order(kgt$n, kgt$lambda), ]
write.csv(kgt, file.path(RES, "sim_main2_kgrid_relerr.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 11. Running time per replicate
## ---------------------------------------------------------------------------
rt <- do.call(rbind, lapply(keys, function(k)
  data.frame(graph = key_g(k), n = key_n(k),
             sec_old = mean(sapply(A[[k]][1:NOLD], `[[`, "seconds")),
             sec_new = mean(sapply(A[[k]][-(1:NOLD)], `[[`, "seconds")))))
write.csv(rt, file.path(RES, "sim_main2_seconds.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(RES, "sessionInfo_main2_summary.txt"))
cat("summary tables written\n")
