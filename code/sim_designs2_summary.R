###############################################################################
## Tables for the further simulation designs (sim_designs2.R)
##
## Author : M. Arashi
## Purpose: reads the replicates in ../results/ext_work2/ and writes
##   sim_designs2_selection.csv  TPR, FDP, MCC, edges and FWER of every
##                               selection rule, with Monte Carlo standard
##                               errors; tie_best marks the rules within two
##                               paired standard errors of the best MCC
##   sim_designs2_ranking.csv    best attainable MCC of every ranking and the
##                               share of true edges among its top |E| pairs
##   sim_designs2_margin.csv     penalties, kappa_hat against kappa_n, the
##                               population margin, the leakage / noise split
##                               and the oracle bias-corrected ranking
##   sim_designs2_sens.csv       the same on the penalty grid of RT-MPS
##   sim_designs2_kfit.csv       is E Q a ridge inverse, and at which kappa
## Usage  : Rscript sim_designs2_summary.R   (from code/)
###############################################################################

OUT  <- ".."
WORK <- file.path(OUT, "results", "ext_work2")
RES  <- file.path(OUT, "results")

files <- list.files(WORK, pattern = "^d[1-4]_.*_rep[0-9]+\\.rds$", full.names = TRUE)
cell_of <- sub("_rep[0-9]+\\.rds$", "", basename(files))

## cell name -> design variables
cell_vars <- function(cell) {
  k <- strsplit(cell, "_")[[1]]
  d <- k[1]
  v <- switch(d,
    d1 = list(graph = k[2], n = as.integer(k[3]), p = 100L, K = NA, dist = "gauss"),
    d2 = list(graph = k[2], n = as.integer(k[3]), p = 100L, K = NA, dist = "gauss"),
    d3 = list(graph = k[2], p = as.integer(k[3]), K = as.integer(k[4]),
              n = as.integer(k[5]), dist = "gauss"),
    d4 = list(graph = k[2], dist = k[3], n = as.integer(k[4]), p = 100L, K = NA))
  data.frame(design = d, graph = v$graph, p = v$p, K = v$K, dist = v$dist,
             n = v$n, stringsAsFactors = FALSE)
}

mse <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

sel_rows <- rank_rows <- marg_rows <- sens_rows <- list()
for (cell in unique(cell_of)) {
  reps <- lapply(files[cell_of == cell], readRDS)
  bad <- sapply(reps, function(r) !is.null(r$error))
  if (any(bad)) cat(cell, ":", sum(bad), "replicates with errors:",
                    unique(sapply(reps[bad], `[[`, "error")), "\n")
  reps <- reps[!bad]
  if (length(reps) < 2) next     # a cell still being computed
  cv <- cell_vars(cell)
  parts <- if (!is.null(reps[[1]]$r)) c(raw = "a", removed = "r") else c(raw = "a")
  for (variant in names(parts)) {
    A <- lapply(reps, `[[`, parts[variant])
    nr <- length(A)
    ## selection rules
    met <- simplify2array(lapply(A, `[[`, "met"))   # method x stat x rep
    mcc <- met[, "MCC", ]
    best <- which.max(rowMeans(mcc))
    tie <- sapply(seq_len(nrow(mcc)), function(i) {
      if (i == best) return(TRUE)
      dd <- mcc[best, ] - mcc[i, ]
      mean(dd) <= 2 * sd(dd) / sqrt(nr)
    })
    fw <- apply(met[, "FDP", , drop = FALSE] > 0, 1, mean)
    sel_rows[[length(sel_rows) + 1]] <- data.frame(cv, variant = variant,
      method = dimnames(met)[[1]], nrep = nr,
      TPR = rowMeans(met[, "TPR", ]), seTPR = apply(met[, "TPR", ], 1, mse),
      FDP = rowMeans(met[, "FDP", ]), seFDP = apply(met[, "FDP", ], 1, mse),
      MCC = rowMeans(mcc), seMCC = apply(mcc, 1, mse),
      NE = rowMeans(met[, "NE", ]), FWER = fw, seFWER = sqrt(fw * (1 - fw) / nr),
      tie_best = tie, row.names = NULL)
    ## rankings
    rk <- simplify2array(lapply(A, `[[`, "rks"))    # method x (best, match) x rep
    bb <- rk[, "best", ]
    bbest <- which.max(rowMeans(bb[rownames(bb) != "Oracle_corr", ]))
    bname <- rownames(bb)[rownames(bb) != "Oracle_corr"][bbest]
    rtie <- sapply(rownames(bb), function(m) {
      dd <- bb[bname, ] - bb[m, ]
      mean(dd) <= 2 * sd(dd) / sqrt(nr)
    })
    rank_rows[[length(rank_rows) + 1]] <- data.frame(cv, variant = variant,
      method = rownames(bb), nrep = nr,
      best_mcc = rowMeans(bb), se_best = apply(bb, 1, mse),
      match_tpr = rowMeans(rk[, "match", ]), se_match = apply(rk[, "match", ], 1, mse),
      tie_best = rtie & rownames(bb) != "Oracle_corr", row.names = NULL)
    ## margin, effective penalty, leakage / noise
    I <- sapply(A, `[[`, "info"); P <- sapply(A, `[[`, "pop")
    rt <- bb["RT", ]; oc <- bb["Oracle_corr", ]
    dk <- I["kappa", ] - I["kappa_mc", ]
    sec <- sapply(reps, `[[`, "seconds")
    mrow <- data.frame(cv, variant = variant, nrep = nr,
      lambda = mean(I["lambda", ]), kappa_hat = mean(I["kappa", ]),
      se_kappa_hat = mse(I["kappa", ]), sd_kappa_hat = sd(I["kappa", ]),
      kappa_mc = mean(I["kappa_mc", ]), kappa_diff = mean(dk), se_kappa_diff = mse(dk),
      kappa_reldiff = mean(dk / I["kappa_mc", ]),
      ratio_kappa_lambda = mean(I["kappa", ] / I["lambda", ]),
      zeta = mean(P["zeta", ]), ell = mean(P["ell", ]),
      margin = mean(P["zeta", ] - P["ell", ]), se_margin = mse(P["zeta", ] - P["ell", ]),
      share_ell_gt_zeta = mean(P["ell", ] > P["zeta", ]),
      pop_mcc = mean(P["pop_mcc", ]), se_pop_mcc = mse(P["pop_mcc", ]),
      rt_mcc = mean(rt), se_rt_mcc = mse(rt),
      leak_error = mean(1 - P["pop_mcc", ]), noise_error = mean(P["pop_mcc", ] - rt),
      oracle_corr_mcc = mean(oc), oracle_gain = mean(oc - rt), se_oracle_gain = mse(oc - rt),
      zeta0 = mean(P["zeta0", ]), ell0 = mean(P["ell0", ]), pop0_mcc = mean(P["pop0_mcc", ]),
      lam_vwp0 = mean(I["lam_vwp0", ]), lam_vwpd = mean(I["lam_vwpd", ]),
      tau_daie = mean(I["tau", ]),
      gn_lambda = mean(I["gn_lambda", ]), kappa_genenet = mean(I["kappa_gn", ]),
      gn_margin = mean(P["gn_zeta", ] - P["gn_ell", ]),
      gn_share_ell_gt_zeta = mean(P["gn_ell", ] > P["gn_zeta", ]),
      gn_pop_mcc = mean(P["gn_pop_mcc", ]),
      beam_alpha = mean(I["beam_alpha", ]), kappa_beam = mean(I["kappa_beam", ]),
      maxdev = mean(I["maxdev", ]), hs_fail = mean(I["hs_ok", ] == 0),
      liu_fail = mean(I["liu_ok", ] == 0), sec_per_rep = mean(sec),
      stringsAsFactors = FALSE)
    if (!is.null(reps[[1]]$model)) {
      M <- sapply(reps, `[[`, "model")
      mrow$cond_med <- median(M["cond", ]); mrow$cond_min <- min(M["cond", ])
      mrow$cond_max <- max(M["cond", ]); mrow$pc_min <- mean(M["pc_min", ])
      mrow$pc_med <- mean(M["pc_med", ]); mrow$pc_max <- mean(M["pc_max", ])
      mrow$n_edges <- mean(M["n_edges", ])
    }
    marg_rows[[length(marg_rows) + 1]] <- mrow
    if (!is.null(A[[1]]$sens)) {
      Sn <- simplify2array(lapply(A, `[[`, "sens"))
      sens_rows[[length(sens_rows) + 1]] <- data.frame(cv, variant = variant,
        mult = c(0.3, 1, 3, 10), kappa = rowMeans(Sn["kappa", , ]),
        rt_mcc = rowMeans(Sn["mcc", , ]), se_rt_mcc = apply(Sn["mcc", , ], 1, mse),
        margin = rowMeans(Sn["zeta", , ] - Sn["ell", , ]),
        share_ell_gt_zeta = rowMeans(Sn["ell", , ] > Sn["zeta", , ]),
        pop_mcc = rowMeans(Sn["pop_mcc", , ]))
    }
  }
}
## rows of different designs can have different columns (d1 has the model
## summaries); missing columns are filled with NA
bind <- function(x) {
  cols <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(d) { for (k in setdiff(cols, names(d))) d[[k]] <- NA; d[cols] })
  y <- do.call(rbind, x); rownames(y) <- NULL; y
}
ord <- function(d) d[order(d$design, d$graph, d$p, d$K, d$dist, d$n), ]
write.csv(ord(bind(sel_rows)), file.path(RES, "sim_designs2_selection.csv"), row.names = FALSE)
write.csv(ord(bind(rank_rows)), file.path(RES, "sim_designs2_ranking.csv"), row.names = FALSE)
write.csv(ord(bind(marg_rows)), file.path(RES, "sim_designs2_margin.csv"), row.names = FALSE)
if (length(sens_rows))
  write.csv(ord(bind(sens_rows)), file.path(RES, "sim_designs2_sens.csv"), row.names = FALSE)

## split of the oracle bias correction (part "oracle" of sim_designs2.R)
of <- list.files(WORK, pattern = "^oracle_d[1-4]_.*\\.rds$", full.names = TRUE)
if (length(of)) {
  ocell <- sub("^oracle_", "", sub("_rep[0-9]+\\.rds$", "", basename(of)))
  orows <- list()
  for (cell in unique(ocell)) {
    O <- lapply(of[ocell == cell], readRDS)
    if (length(O) < 2) next
    B <- sapply(O, function(o) o$rks[, "best"])
    chk <- max(abs(sapply(O, `[[`, "check")))
    g <- function(k) B[k, ] - B["RT", ]
    orows[[length(orows) + 1]] <- data.frame(cell_vars(cell), nrep = ncol(B),
      rt_mcc = mean(B["RT", ]),
      gain_leak_only = mean(g("leak_only")), se_leak_only = mse(g("leak_only")),
      gain_shrink_only = mean(g("shrink_only")), se_shrink_only = mse(g("shrink_only")),
      gain_full = mean(g("full")), se_full = mse(g("full")),
      max_check = chk)
  }
  write.csv(ord(bind(orows)), file.path(RES, "sim_designs2_oracle.csv"), row.names = FALSE)
}

kf <- list.files(WORK, pattern = "^kfit_.*\\.rds$", full.names = TRUE)
if (length(kf)) {
  K <- bind(lapply(kf, readRDS))
  write.csv(K[order(K$design, K$graph, K$p, K$K, K$dist, K$scores, K$n, K$lambda), ],
            file.path(RES, "sim_designs2_kfit.csv"), row.names = FALSE)
}
writeLines(capture.output(sessionInfo()),
           file.path(RES, "sessionInfo_designs2_summary.txt"))
cat("tables written\n")
