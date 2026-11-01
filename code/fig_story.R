###############################################################################
## Simulation figures of Section 5 (Figures 2 to 4) of
##   "What does a ridge network estimate? Effective penalty and leakage in
##    ridge precision matrix estimation"
##
## Author : M. Arashi
## Purpose: three figures that summarize the simulation tables; the exact
##          values are in the supplement.
##          Figure 2 (fig_sim_loss.pdf): where the accuracy of the sample
##              ridge ranking goes. For each setting, 1 - best attainable MCC
##              of RT is split into the gain from removing the population
##              leakage from the non-edges (leakage) and the loss that is left
##              after that correction (noise).
##          Figure 3 (fig_sim_mcc.pdf): MCC of the selected graphs and best
##              attainable MCC of the rankings, by method, graph and n.
##          Figure 4 (fig_sim_designs.pdf): designs D1 to D4, relative error
##              of the estimated effective penalty and best attainable MCC of
##              RT against the MB path.
##          Methods carry one colour per family in Figures 3 and 4: ridge
##          rules blue, ridge-type rules calibrated by false discovery rates
##          vermillion, lasso and nodewise methods green, marginal grey.
## Input  : ../results/sim_main2_leakage_oracle.csv,
##          ../results/sim_main2_selection.csv, ../results/sim_main2_ranking.csv,
##          ../results/sim_designs2_margin.csv, ../results/sim_designs2_ranking.csv,
##          ../results/sim_designs2_paths.csv,
##          ../results/sim_designs2_factor_fair.csv
## Output : ../figs/fig_sim_loss.pdf, fig_sim_mcc.pdf, fig_sim_designs.pdf,
##          copied to ../figures/ by make_figures.R
## Usage  : Rscript fig_story.R   (run from code/); also called from
##          make_figures.R
###############################################################################

OUT <- ".."
RES <- file.path(OUT, "results")
FIG <- file.path(OUT, "figs")
dir.create(FIG, showWarnings = FALSE)

## printed size: text width of the journal is about 5.6 in, text at 9 pt
open_fig <- function(name, width, height)
  cairo_pdf(file.path(FIG, name), width = width, height = height,
            pointsize = 9, family = "sans")

okabe <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
           blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7",
           black = "#000000")
fam_col <- c(ridge = okabe[["blue"]], fdr = okabe[["vermillion"]],
             lasso = okabe[["green"]], marg = "gray35")
grid_col <- "gray88"
graphs <- c("band", "random", "hub", "scale-free")
gl <- c(band = "Band", random = "Random", hub = "Hub", "scale-free" = "Scale-free")
n_grid <- c(50, 100, 200, 400)
pch_n <- c(21, 22, 24, 23)                       # n = 50, 100, 200, 400


## horizontal key drawn by hand: symbols and labels in one row, starting at x0
hkey <- function(x0, y, labels, pch = NA, bg = NA, col = "gray20", fill = NA,
                 cex = 0.9, gap = 0.025) {
  x <- x0
  for (k in seq_along(labels)) {
    if (!is.na(fill[k])) {
      rect(x, y - 0.18, x + 0.018, y + 0.18, col = fill[k], border = NA)
    } else {
      points(x + 0.009, y, pch = pch[k], bg = bg[k], col = col, cex = 1)
    }
    text(x + 0.026, y, labels[[k]], adj = 0, cex = cex)
    x <- x + 0.026 + strwidth(labels[[k]], cex = cex) + gap
  }
  invisible(x)
}

## colour of the same hue, from light (small n) to full (large n)
shade <- function(col, k = 4) colorRampPalette(c("white", col))(k + 2)[-(1:2)]

## ===========================================================================
## Figure 2: where the accuracy goes
## ===========================================================================
lo <- read.csv(file.path(RES, "sim_main2_leakage_oracle.csv"))
lo$leak <- lo$ORC_leak - lo$RT                   # = d_leak
lo$noise <- 1 - lo$ORC_leak
stopifnot(abs(lo$leak - lo$d_leak) < 1e-12)

## row positions: groups of bars with a gap between groups
bar_rows <- function(groups, n_by_group) {
  y <- numeric(0); g <- character(0); n <- numeric(0); head <- numeric(0)
  pos <- 0
  for (k in seq_along(groups)) {
    pos <- pos + 1.0                             # room for the group label
    head <- c(head, pos)
    for (nn in n_by_group[[k]]) {
      pos <- pos + 1
      y <- c(y, pos); g <- c(g, groups[k]); n <- c(n, nn)
    }
    pos <- pos + 0.35
  }
  list(y = y, g = g, n = n, head = head, max = pos)
}

draw_loss <- function(d, rows, ymax, group_lab) {
  plot(NA, xlim = c(0, 1), ylim = c(ymax, 0.4), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  abline(v = seq(0, 1, 0.2), col = grid_col, lwd = 0.6)
  for (i in seq_along(rows$y)) {
    r <- d[d$graph == rows$g[i] & d$n == rows$n[i], ]
    stopifnot(nrow(r) == 1)
    y <- rows$y[i]; h <- 0.36
    rect(0, y - h, r$leak, y + h, col = okabe[["purple"]], border = NA)
    rect(r$leak, y - h, r$leak + r$noise, y + h, col = "gray78", border = NA)
    ## long bars carry their value inside the bar, so it cannot touch the column
    tot <- r$leak + r$noise
    if (tot > 0.8) {
      text(tot - 0.015, y, sprintf("%.2f", 1 - r$RT), adj = 1, cex = 0.9, col = "gray20")
    } else {
      text(tot + 0.015, y, sprintf("%.2f", 1 - r$RT), adj = 0, cex = 0.9, col = "gray25")
    }
    ## the leakage part as a number, in a column right of the panel
    text(1.2, y, sprintf("%.3f", r$leak), adj = 1, xpd = NA, cex = 0.9,
         col = "#9E3D7A")
  }
  text(1.2, rows$head[1], "leak.", adj = 1, xpd = NA, cex = 0.9,
       col = "#9E3D7A", font = 2)
  axis(2, at = rows$y, labels = rows$n, las = 1, tick = FALSE, line = -0.35,
       cex.axis = 0.9)
  text(-0.01, rows$head, group_lab, adj = c(1, 0.5), xpd = NA, font = 2,
       cex = 1)
  axis(1, at = seq(0, 1, 0.2), labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"),
       mgp = c(2, 0.45, 0), tcl = -0.25, cex.axis = 0.9)
}

main_rows <- bar_rows(graphs, rep(list(n_grid), 4))
hd_groups <- c("band, p = 400", "hub, p = 400", "random, p = 400")
hd_rows <- bar_rows(hd_groups, rep(list(c(50, 100)), 3))
ymax <- max(main_rows$max, hd_rows$max)

open_fig("fig_sim_loss.pdf", 5.6, 3.5)
layout(matrix(1:2, 1), widths = c(1, 1))
par(cex = 1)
par(mar = c(2.9, 4.4, 1.9, 3.0), mgp = c(1.8, 0.5, 0), las = 1)
draw_loss(lo, main_rows, ymax, gl[graphs])
mtext(expression("(a) main design, " * italic(p) == 100), side = 3, line = 0.55,
      adj = 0, at = -0.3, cex = 1)
draw_loss(lo, hd_rows, ymax, c("Band", "Hub", "Random"))
mtext(expression("(b) block with " * italic(p) == 400), side = 3, line = 0.55,
      adj = 0, at = -0.3, cex = 1)
## one x label for both panels, and the key in the free space of panel (b)
mtext("1 − best attainable MCC of the sample ridge ranking RT", side = 1,
      line = 1.75, at = -0.25, cex = 0.9)
key_y <- hd_rows$max + c(2.2, 4.9)
key_txt <- list(c("leakage: gain from removing",
                  "the population leakage"),
                c("noise: loss left after", "removing the leakage"))
key_col <- c(okabe[["purple"]], "gray78")
for (k in 1:2) {
  rect(0.02, key_y[k] - 0.36, 0.09, key_y[k] + 0.36, col = key_col[k], border = NA)
  text(0.12, key_y[k], key_txt[[k]][1], adj = 0, cex = 0.9)
  text(0.12, key_y[k] + 0.95, key_txt[[k]][2], adj = 0, cex = 0.9)
}
dev.off()

## ===========================================================================
## Figure 3: who wins
## ===========================================================================
se <- read.csv(file.path(RES, "sim_main2_selection.csv"))
se <- se[se$reps == 500, ]
rk <- read.csv(file.path(RES, "sim_main2_ranking.csv"))
rk <- rk[rk$reps == 500, ]

## rows from top to bottom: label, column in the csv, family
sel_methods <- data.frame(
  lab = c("MB", "GLASSO", "Liu", "HS", "beam", "GeneNet", "vWP0-Fdr",
          "RS", "vWPd", "vWP0", "RT"),
  key = c("MB", "GLASSO", "Liu2013", "HS_BH10", "beam_BH10", "GeneNet_lFDR80",
          "vWP_Fdr10", "RS_lasso", "vWPt_lFDR90", "vWP_lFDR", "RT_lFDR"),
  fam = c("lasso", "lasso", "lasso", "fdr", "fdr", "fdr", "fdr",
          "ridge", "ridge", "ridge", "ridge"), stringsAsFactors = FALSE)
rank_methods <- data.frame(
  lab = c("MB", "GLASSO", "Marginal", "beam", "GeneNet", "vWP0", "RT, stab.", "RT"),
  key = c("MB", "GLASSO", "Marginal", "beam", "GeneNet", "vWP", "RT_stab", "RT"),
  fam = c("lasso", "lasso", "marg", "fdr", "fdr", "ridge", "ridge", "ridge"),
  stringsAsFactors = FALSE)

get_val <- function(tab, valcol, key, g, n) {
  v <- tab[tab$method == key & tab$graph == g & tab$n == n, valcol]
  stopifnot(length(v) == 1)
  v
}

draw_dots <- function(tab, valcol, meth, g, first_col, last_row) {
  m <- nrow(meth)
  plot(NA, xlim = c(-0.03, 1.03), ylim = c(m + 0.6, 0.4), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  abline(v = seq(0, 1, 0.25), col = grid_col, lwd = 0.6)
  ## thin separators between the method families
  brk <- which(meth$fam[-1] != meth$fam[-m])
  abline(h = brk + 0.5, col = "gray70", lwd = 0.6)
  for (i in seq_len(m)) {
    v <- sapply(n_grid, function(n) get_val(tab, valcol, meth$key[i], g, n))
    segments(min(v), i, max(v), i, col = "gray60", lwd = 0.8)
    cols <- shade(fam_col[[meth$fam[i]]])
    points(v, rep(i, 4), pch = pch_n, bg = cols, col = fam_col[[meth$fam[i]]],
           cex = 0.95, lwd = 0.7)
  }
  box(col = "gray55", lwd = 0.6)
  if (first_col)
    axis(2, at = seq_len(m), labels = meth$lab, las = 1, tick = FALSE,
         line = -0.65, cex.axis = 0.9)
  axis(1, at = seq(0, 1, 0.25), labels = if (last_row)
         c("0", "", "0.5", "", "1") else FALSE,
       mgp = c(2, 0.4, 0), tcl = -0.2, cex.axis = 0.9)
}

open_fig("fig_sim_mcc.pdf", 5.6, 4.75)
layout(matrix(c(1, 1, 1, 1, 2:9), 3, byrow = TRUE),
       widths = c(1.88, 1.24, 1.24, 1.24),
       heights = c(1.7, nrow(sel_methods) + 1.6, nrow(rank_methods) + 2.3))
par(cex = 1)                                     # layout() would shrink the text
## key: shapes for n, colours for the method families
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(c(0, 1), c(-1, 1))
hkey(0.02, 0.45, list(expression(italic(n) == 50), "100", "200", "400"),
     pch = pch_n, bg = shade("gray30"))
hkey(0.02, -0.5, c("ridge rules", "FDR-calibrated ridge type", "lasso type",
                   "marginal"),
     fill = fam_col[c("ridge", "fdr", "lasso", "marg")])
for (j in seq_along(graphs)) {
  par(mar = c(1.3, if (j == 1) 4.6 else 0.35, 1.5, 0.35), mgp = c(1.6, 0.4, 0))
  draw_dots(se, "MCC", sel_methods, graphs[j], j == 1, FALSE)
  mtext(gl[[graphs[j]]], side = 3, line = 0.25, cex = 1, font = 2)
  if (j == 1) mtext("(a) selection", side = 3, line = 0.25,
                    at = grconvertX(0.01, "ndc", "user"),
                    adj = 0, cex = 0.9)
}
for (j in seq_along(graphs)) {
  par(mar = c(2.6, if (j == 1) 4.6 else 0.35, 1.5, 0.35), mgp = c(1.6, 0.4, 0))
  draw_dots(rk, "best", rank_methods, graphs[j], j == 1, TRUE)
  if (j == 1) mtext("(b) ranking", side = 3, line = 0.25,
                    at = grconvertX(0.01, "ndc", "user"),
                    adj = 0, cex = 0.9)
  if (j == 2) mtext("MCC of the selected graph (a); best attainable MCC (b)",
                    side = 1, line = 1.5, at = 1.05, cex = 0.9)
}
dev.off()

## ===========================================================================
## Figure 4: robustness across designs
## ===========================================================================
dm <- read.csv(file.path(RES, "sim_designs2_margin.csv"))
dm <- dm[dm$variant == "raw", ]
dr <- read.csv(file.path(RES, "sim_designs2_ranking.csv"))
dp <- read.csv(file.path(RES, "sim_designs2_paths.csv"))
ff <- read.csv(file.path(RES, "sim_designs2_factor_fair.csv"))

pch_g <- c(band = 21, random = 24, hub = 22, "scale-free" = 23)
col_g <- "gray20"
designs <- c("d1", "d2", "d3", "d4")
dlab <- list(d1 = "D1 heterogeneous weights", d2 = expression("D2 large " * italic(n)),
             d3 = "D3 factors", d4 = "D4 heavy tails")
cell_id <- function(d) paste(d$design, d$graph, d$p, d$K, d$dist, d$n)

## best attainable MCC of RT (rankings) and of the MB path (rescored path)
rt <- dr[dr$method == "RT", ]
mb <- dp[dp$method == "MB", ]
pair <- merge(rt[, c("design", "graph", "p", "K", "dist", "n", "variant", "best_mcc")],
              mb[, c("design", "graph", "p", "K", "dist", "n", "variant", "best_graph")],
              by = c("design", "graph", "p", "K", "dist", "n", "variant"))
## D3 after removing the components: the factor-adjusted lasso comparison
ffr <- ff[ff$variant == "removed" & ff$method %in% c("RT_rank", "MB_path"), ]
chk <- merge(pair[pair$design == "d3" & pair$variant == "removed", ],
             reshape(ffr[, c("graph", "p", "K", "n", "method", "MCC")],
                     idvar = c("graph", "p", "K", "n"), timevar = "method",
                     direction = "wide"), by = c("graph", "p", "K", "n"))
stopifnot(nrow(chk) == 24, abs(chk$best_mcc - chk$MCC.RT_rank) < 1e-9,
          abs(chk$best_graph - chk$MCC.MB_path) < 1e-9)

open_fig("fig_sim_designs.pdf", 5.6, 3.7)
layout(matrix(c(1, 1, 1, 1, 2:9), 3, byrow = TRUE), heights = c(0.3, 1.0, 1.25),
       widths = c(1.805, 1.265, 1.265, 1.265))
par(cex = 1)                                     # layout() would shrink the text
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(c(0, 1), c(-1, 1))
hkey(0.02, 0, as.list(unname(gl)), pch = pch_g[graphs], bg = rep("gray55", 4))
hkey(0.5, 0, list("D3 raw data", "D3 top components removed"), pch = c(21, 21),
     bg = c("white", "gray55"))

## top row: relative error of kappa-hat (%) against n
for (j in seq_along(designs)) {
  d <- dm[dm$design == designs[j], ]
  par(mar = c(2.3, if (j == 1) 4.2 else 0.6, 1.5, 0.6), mgp = c(1.8, 0.4, 0),
      las = 1)
  xl <- range(d$n) * c(0.8, 1.25)
  plot(NA, xlim = xl, ylim = c(-8.5, 1), log = "x", axes = FALSE,
       xlab = "", ylab = "")
  abline(h = c(-8, -6, -4, -2), col = grid_col, lwd = 0.6)
  abline(h = 0, col = "gray55", lwd = 0.8)
  for (cfg in unique(paste(d$graph, d$p, d$K, d$dist))) {
    e <- d[paste(d$graph, d$p, d$K, d$dist) == cfg, ]
    e <- e[order(e$n), ]
    lines(e$n, 100 * e$kappa_reldiff, col = "gray65", lwd = 0.7)
  }
  points(d$n, 100 * d$kappa_reldiff, pch = pch_g[d$graph], bg = "gray55",
         col = col_g, cex = 0.9, lwd = 0.6)
  box(col = "gray55", lwd = 0.6)
  nn <- sort(unique(d$n))
  if (length(nn) == 5) nn_lab <- ifelse(nn %in% c(50, 200, 800), nn, NA) else nn_lab <- nn
  axis(1, at = nn, labels = nn_lab, mgp = c(2, 0.35, 0), tcl = -0.2,
       cex.axis = 0.9, gap.axis = 0.2)
  if (j == 1) axis(2, at = c(-8, -6, -4, -2, 0), tcl = -0.2, cex.axis = 0.9)
  if (j == 1) mtext("relative error (%)", side = 2,
                    line = 2.4, las = 0, cex = 0.9)
  mtext(dlab[[designs[j]]], side = 3, line = 0.25, cex = 0.9)
  mtext(expression(italic(n)), side = 1, line = 1.2, cex = 0.9)
}

## bottom row: best MCC of RT (x) against the MB path (y)
for (j in seq_along(designs)) {
  d <- pair[pair$design == designs[j], ]
  par(mar = c(2.6, if (j == 1) 4.2 else 0.6, 0.6, 0.6), mgp = c(1.8, 0.4, 0),
      las = 1, pty = "m")
  plot(NA, xlim = c(0, 1.02), ylim = c(0, 1.02), xaxs = "i", yaxs = "i",
       axes = FALSE, xlab = "", ylab = "")
  polygon(c(0, 0, 1.02), c(0, 1.02, 1.02), col = adjustcolor(fam_col[["lasso"]], 0.10),
          border = NA)
  polygon(c(0, 1.02, 1.02), c(0, 0, 1.02), col = adjustcolor(fam_col[["ridge"]], 0.10),
          border = NA)
  abline(v = seq(0.25, 1, 0.25), h = seq(0.25, 1, 0.25), col = "white", lwd = 0.6)
  abline(0, 1, col = "gray45", lwd = 0.8)
  if (designs[j] == "d2") {                     # D2 has one point: room for the key
    text(0.05, 0.9, "MB path
higher", adj = c(0, 1), cex = 0.9, col = "#00664A")
    text(0.97, 0.08, "RT higher", adj = c(1, 0), cex = 0.9, col = "#004C77")
  }
  fillc <- ifelse(d$variant == "removed" | d$design != "d3", "gray55", "white")
  points(d$best_mcc, d$best_graph, pch = pch_g[d$graph], bg = fillc, col = col_g,
         cex = 0.9, lwd = 0.6)
  box(col = "gray55", lwd = 0.6)
  axis(1, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), mgp = c(2, 0.35, 0),
       tcl = -0.2, cex.axis = 0.9)
  if (j == 1) axis(2, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), tcl = -0.2,
                   cex.axis = 0.9)
  if (j == 1) mtext("best MCC, MB path", side = 2, line = 2.4, las = 0, cex = 0.9)
  mtext("best MCC, RT", side = 1, line = 1.3, cex = 0.9)
}
dev.off()

cat("story figures written: fig_sim_loss.pdf fig_sim_mcc.pdf fig_sim_designs.pdf\n")
