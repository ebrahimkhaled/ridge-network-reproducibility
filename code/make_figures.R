###############################################################################
## Figures for
##   "What does a ridge network estimate? Effective penalty, leakage and
##    support recovery for ridge precision matrix estimators"
##
## Author : M. Arashi
## Purpose: draws all figures from the saved results, so they can be redrawn
##          without rerunning any computation. cairo_pdf() is used, so the
##          fonts are embedded and Greek letters are kept. Lines for different
##          n or lambda carry both a colour (viridis, colour-blind safe) and a
##          point shape; categories use the Okabe-Ito palette.
## Input  : ../results/*.csv, ../results/data_selections.rds, the stock data
##          of the huge package (sector labels)
## Output : ../figs/*.pdf, copied to ../figures/ (the folder the manuscript
##          reads)
## Usage  : Rscript make_figures.R   (run from code/)
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages({ library(igraph); library(viridisLite) })
OUT <- ".."
RES <- file.path(OUT, "results")
FIG <- file.path(OUT, "figs")
dir.create(FIG, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)

okabe <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
           yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00",
           purple = "#CC79A7", black = "#000000")
seq4 <- viridis(4, begin = 0, end = 0.8)        # n or lambda, light = small
pch4 <- c(16, 15, 17, 18)
graphs <- c("band", "random", "hub", "scale-free")
gl <- c(band = "Band", random = "Random", hub = "Hub", "scale-free" = "Scale-free")
n_grid <- c(50, 100, 200, 400)
open_pdf <- function(name, width, height)
  cairo_pdf(file.path(FIG, name), width = width, height = height, pointsize = 10,
            family = "sans")
## labels 10^k on a log axis
pow10_axis <- function(side, lim) {
  k <- floor(log10(lim[1])):ceiling(log10(lim[2]))
  axis(side, at = 10^k, labels = as.expression(lapply(k, function(e) bquote(10^.(e)))),
       las = 1)
}

## ---------------------------------------------------------------------------
## Figure 1: effective penalty, band model, p = 200
## ---------------------------------------------------------------------------
kg <- read.csv(file.path(RES, "sim_kappa_grid.csv"))
open_pdf("fig_sim_kappa.pdf", 6.5, 3.6)
par(mar = c(3.8, 4.2, 3.4, 0.8), mgp = c(2.4, 0.7, 0), las = 1)
lams <- sort(unique(kg$lambda))
r <- 200 / kg$n
plot(NA, xlim = range(r), ylim = c(0, max(kg$kappa_mc, kg$kappa_hat) * 1.05),
     log = "x", xlab = expression(italic(p) / italic(n)),
     ylab = "effective penalty", xaxt = "n")
axis(1, at = sort(unique(r)), labels = c("1/4", "1/2", "1", "2", "4"))
for (k in seq_along(lams)) {
  d <- kg[kg$lambda == lams[k], ]; d <- d[order(200 / d$n), ]
  lines(200 / d$n, d$kappa_mc, col = seq4[k], lwd = 1.6)
  points(200 / d$n, d$kappa_hat, col = seq4[k], pch = pch4[k], cex = 1.1)
  lines(200 / d$n, d$floor, col = seq4[k], lty = 3, lwd = 1.3)
}
## two legends above the plot: nominal penalty (colour and shape) and the
## meaning of lines and points
legend("bottomleft", inset = c(0, 1.02), xpd = NA, horiz = TRUE, bty = "n",
       legend = as.expression(lapply(lams, function(l) bquote(lambda == .(l)))),
       col = seq4, pch = pch4, lwd = 1.6, cex = 0.9)
legend("bottomleft", inset = c(0, 1.1), xpd = NA, horiz = TRUE, bty = "n",
       legend = c(expression("Monte Carlo " * kappa[italic(n)]),
                  expression("mean of " * hat(kappa)[lambda]),
                  "lower bound"),
       lty = c(1, NA, 3), pch = c(NA, 16, NA), lwd = c(1.6, NA, 1.3),
       col = "gray25", cex = 0.9)
dev.off()

## ---------------------------------------------------------------------------
## Figure 2 (as published): mean |partial correlation| by graph distance
## ---------------------------------------------------------------------------
lk <- read.csv(file.path(RES, "sim_leakage.csv"))
open_pdf("fig_sim_leakage.pdf", 6.5, 5.4)
par(mfrow = c(2, 2), mar = c(3.6, 4.2, 1.6, 0.6), mgp = c(2.4, 0.7, 0),
    oma = c(0, 0, 3.6, 0), las = 1)
for (g in graphs) {
  d <- lk[lk$graph == g, ]
  plot(NA, xlim = c(1, 4), ylim = c(0, max(d$observed, d$predicted, na.rm = TRUE)),
       xaxt = "n", xlab = "graph distance", ylab = "mean |partial correlation|",
       main = gl[g], font.main = 1)
  axis(1, at = 1:4, labels = c("1", "2", "3", "4+"))
  for (k in seq_along(n_grid)) {
    e <- d[d$n == n_grid[k], ]
    ok <- !is.na(e$observed)
    lines((1:4)[ok], e$observed[ok], col = seq4[k], lwd = 1.5)
    points((1:4)[ok], e$observed[ok], col = seq4[k], pch = pch4[k])
    lines((1:4)[ok], e$predicted[ok], col = seq4[k], lwd = 1.5, lty = 2)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("top", horiz = TRUE, bty = "n", legend = paste0("n = ", n_grid),
       col = seq4, pch = pch4, lwd = 1.5, cex = 0.9)
legend("top", inset = c(0, 0.045), horiz = TRUE, bty = "n",
       legend = c("observed", "population ridge network"), lty = c(1, 2),
       lwd = 1.5, col = "gray25", cex = 0.9)
dev.off()

## ---------------------------------------------------------------------------
## New Figure 2: mean signed partial correlation of every pair against the
## population ridge network and against the partial correlation of Theta
## ---------------------------------------------------------------------------
pr <- read.csv(file.path(RES, "sim_pairs_signed.csv"))
pr$dclass <- pmin(pr$distance, 4)
dcol <- c(okabe["black"], okabe["vermillion"], okabe["blue"], okabe["skyblue"])
dpch <- c(16, 17, 15, 1)
signed_panels <- function(nsel, file) {
  open_pdf(file, 6.5, 3.9)
  par(mfrow = c(2, 4), mar = c(3.4, 4.3, 1.6, 0.4), mgp = c(2.9, 0.6, 0),
      oma = c(0, 0, 1.8, 0), las = 1, cex = 0.75)
  for (xvar in c("rho_kappa_mean", "rho_theta")) for (g in graphs) {
    d <- pr[pr$graph == g & pr$n == nsel, ]
    d <- d[order(-d$dclass), ]                   # far pairs drawn first
    lim <- range(d$rho_hat_mean, d[[xvar]])
    plot(d[[xvar]], d$rho_hat_mean, col = dcol[d$dclass], pch = dpch[d$dclass],
         cex = 0.5, xlim = lim, ylim = lim, main = gl[g], font.main = 1,
         xlab = if (xvar == "rho_theta") expression("partial correlation of " * Theta)
                else expression(rho^(kappa[italic(n)])),
         ylab = expression("mean of " * hat(rho)))
    abline(0, 1, col = "gray60")
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot.new()
  legend("top", horiz = TRUE, bty = "n", cex = 0.85,
         legend = c("edges", "distance 2", "distance 3", "distance 4+"),
         col = dcol, pch = dpch)
  dev.off()
}
signed_panels(100, "fig_sim_leakage_signed.pdf")
for (nn in n_grid) signed_panels(nn, sprintf("fig_sim_leakage_signed_n%d.pdf", nn))

## ---------------------------------------------------------------------------
## Figure 3: penalty sensitivity of the best attainable MCC
## ---------------------------------------------------------------------------
sens <- read.csv(file.path(RES, "sim_sensitivity.csv"))
open_pdf("fig_sim_sensitivity.pdf", 6.5, 5.2)
par(mfrow = c(2, 2), mar = c(3.6, 4.2, 1.6, 0.6), mgp = c(2.4, 0.7, 0),
    oma = c(0, 0, 1.4, 0), las = 1)
for (g in graphs) {
  d <- sens[sens$graph == g, ]
  plot(NA, xlim = range(sens$kappa), ylim = c(0, 1), log = "x", xaxt = "n",
       xlab = "mean effective penalty", ylab = "best attainable MCC",
       main = gl[g], font.main = 1)
  axis(1, at = c(0.1, 0.3, 1, 3, 10, 30), labels = c("0.1", "0.3", "1", "3", "10", "30"))
  for (k in seq_along(n_grid)) {
    e <- d[d$n == n_grid[k], ]
    lines(e$kappa, e$mcc, col = seq4[k], lwd = 1.5)
    points(e$kappa, e$mcc, col = seq4[k], pch = pch4[k])
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("top", horiz = TRUE, bty = "n", legend = paste0("n = ", n_grid),
       col = seq4, pch = pch4, lwd = 1.5, cex = 0.9)
dev.off()

## ---------------------------------------------------------------------------
## Figure 4: effective penalty in the stock data
## ---------------------------------------------------------------------------
curves <- read.csv(file.path(RES, "data_kappa_curve.csv"))
marg <- read.csv(file.path(RES, "data_marginal.csv"))
wname <- c(paste0("W", 1:5), "Full")
wcol <- c(viridis(5, end = 0.85), "black")
wpch <- c(16, 15, 17, 18, 25, NA)
open_pdf("fig_data_kappa.pdf", 6.5, 3.8)
par(mfrow = c(1, 2), mar = c(4, 4.6, 1, 0.6), mgp = c(2.9, 0.7, 0),
    oma = c(0, 0, 2.2, 0), las = 1)
xl <- range(curves$lambda); yl <- c(1e-3, 200)
plot(NA, xlim = xl, ylim = yl, log = "xy", xaxt = "n", yaxt = "n",
     xlab = expression("nominal penalty " * lambda),
     ylab = expression("effective penalty " * hat(kappa)[lambda]))
pow10_axis(1, xl); pow10_axis(2, yl)
abline(0, 1, col = "gray75")
for (w in seq_along(wname)) {
  cw <- curves[curves$window == wname[w], ]
  lines(cw$lambda, cw$kappa, col = wcol[w], lty = if (w == 6) 2 else 1, lwd = 1.4)
  if (w < 6) { s <- seq(1, nrow(cw), by = 12)
    points(cw$lambda[s], cw$kappa[s], col = wcol[w], pch = wpch[w], bg = wcol[w], cex = 0.7) }
}
xl2 <- range(marg$kappa)
plot(NA, xlim = xl2, ylim = c(0, 1), log = "x", xaxt = "n",
     xlab = expression("effective penalty " * hat(kappa)[lambda]),
     ylab = "Spearman correlation")
pow10_axis(1, xl2)
for (w in seq_along(wname)) {
  mw <- marg[marg$window == wname[w], ]
  lines(mw$kappa, mw$spearman, col = wcol[w], lty = if (w == 6) 2 else 1, lwd = 1.4)
  if (w < 6) { s <- seq(1, nrow(mw), by = 4)
    points(mw$kappa[s], mw$spearman[s], col = wcol[w], pch = wpch[w], bg = wcol[w], cex = 0.7) }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("top", horiz = TRUE, bty = "n", cex = 0.85,
       legend = c(paste("window", 1:5), "full sample"),
       col = wcol, pch = wpch, pt.bg = wcol, lty = c(1, 1, 1, 1, 1, 2), lwd = 1.4)
legend("top", inset = c(0, 0.06), horiz = TRUE, bty = "n", cex = 0.85,
       legend = c("yearly windows: p/n = 1.8", "full sample: p/n = 0.36",
                  "identity (left)"),
       col = c(NA, NA, "gray75"), lty = c(NA, NA, 1), lwd = 1.4)
dev.off()

## ---------------------------------------------------------------------------
## Figure 3 of the shortened paper: left panel of Figure 4 alone (effective
## penalty against nominal penalty); the two-panel version is Figure S4 of
## the supplement
## ---------------------------------------------------------------------------
open_pdf("fig_data_kappa_left.pdf", 4.6, 3.9)
par(mar = c(4, 4.6, 3.6, 0.8), mgp = c(2.9, 0.7, 0), las = 1)
plot(NA, xlim = xl, ylim = yl, log = "xy", xaxt = "n", yaxt = "n",
     xlab = expression("nominal penalty " * lambda),
     ylab = expression("effective penalty " * hat(kappa)[lambda]))
pow10_axis(1, xl); pow10_axis(2, yl)
abline(0, 1, col = "gray75")
for (w in seq_along(wname)) {
  cw <- curves[curves$window == wname[w], ]
  lines(cw$lambda, cw$kappa, col = wcol[w], lty = if (w == 6) 2 else 1, lwd = 1.4)
  if (w < 6) { s <- seq(1, nrow(cw), by = 12)
    points(cw$lambda[s], cw$kappa[s], col = wcol[w], pch = wpch[w], bg = wcol[w], cex = 0.7) }
}
legend("bottomleft", inset = c(0, 1.01), xpd = NA, ncol = 3, bty = "n", cex = 0.8,
       legend = c(paste("window", 1:5), "full sample"),
       col = wcol, pch = wpch, pt.bg = wcol, lty = c(1, 1, 1, 1, 1, 2), lwd = 1.4)
dev.off()

## ---------------------------------------------------------------------------
## Figure 5: sector concordance and stability (horizontal bars; ridge-based
## and lasso-based methods in different colours)
## ---------------------------------------------------------------------------
comp <- read.csv(file.path(RES, "data_methods.csv"))
lab <- c(vWP_lFDR = "vWP-lFDR", RT_lFDR = "RT-lFDR", RE_lFDR = "RE-lFDR",
         RS_lasso = "RS", RT_MPS = "RT-MPS", GLASSO = "GLASSO", MB = "MB")
comp <- comp[match(names(lab), comp$method), ]
bcol <- ifelse(comp$method %in% c("GLASSO", "MB"), okabe["vermillion"], okabe["blue"])
open_pdf("fig_data_sector.pdf", 6.5, 3.3)
par(mfrow = c(1, 2), mar = c(4, 5.6, 1, 0.8), mgp = c(2.4, 0.7, 0),
    oma = c(0, 0, 1.6, 0), las = 1)
ord <- rev(seq_len(nrow(comp)))                  # first method at the top
barplot(comp$within_sector[ord], names.arg = lab[ord], horiz = TRUE,
        xlim = c(0, 1), col = bcol[ord], border = NA,
        xlab = "within-sector share of edges")
abline(v = comp$baseline[1], lty = 2)
barplot(comp$jaccard_consecutive[ord], names.arg = lab[ord], horiz = TRUE,
        xlim = c(0, 0.6), col = bcol[ord], border = NA,
        xlab = "Jaccard index, consecutive windows")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("top", horiz = TRUE, bty = "n", cex = 0.85,
       legend = c("ridge-based", "lasso-based", "share of within-sector pairs"),
       fill = c(okabe["blue"], okabe["vermillion"], NA),
       border = c(NA, NA, NA), lty = c(NA, NA, 2), x.intersp = 0.8)
dev.off()

## ---------------------------------------------------------------------------
## Figure 6: RT-MPS network of the first window, coloured by GICS sector
## ---------------------------------------------------------------------------
data(stockdata, package = "huge", envir = environment())
sector <- stockdata$info[, 2]
sel_all <- readRDS(file.path(RES, "data_selections.rds"))
A1 <- sel_all[["W1"]][["RT_MPS"]] * 1
keep <- which(rowSums(A1) > 0)
g1 <- graph_from_adjacency_matrix(A1[keep, keep], mode = "undirected")
sec_levels <- sort(unique(sector))
## eight Okabe-Ito colours for the eight sectors with most plotted stocks,
## grey and white (black outline) for the remaining two
cnt <- table(factor(sector[keep], levels = sec_levels))
rank_sec <- names(sort(cnt, decreasing = TRUE))
pal_by_rank <- c(okabe[c("blue", "orange", "green", "vermillion", "skyblue",
                         "purple", "yellow", "black")], "#999999", "#FFFFFF")
pal <- setNames(pal_by_rank, rank_sec)[sec_levels]
vcol <- pal[sector[keep]]
vframe <- ifelse(vcol == "#FFFFFF", "black", NA)
set.seed(1)
lay <- layout_with_fr(g1)
open_pdf("fig_data_network.pdf", 6.5, 4.8)
layout(matrix(1:2, 1), widths = c(4.6, 1.9))
par(mar = c(0, 0, 0, 0))
plot(g1, vertex.size = 3, vertex.label = NA, vertex.color = vcol,
     vertex.frame.color = vframe, edge.color = "gray65", edge.width = 0.6,
     layout = lay)
plot.new()
legend("left", legend = paste0(sec_levels, " (", cnt[sec_levels], ")"),
       pch = 21, pt.bg = pal, col = ifelse(pal == "#FFFFFF", "black", pal),
       pt.cex = 1.3, bty = "n", cex = 0.75, title = "GICS sector (stocks)",
       title.adj = 0)
dev.off()

## copy to the folder the manuscript reads
figs <- list.files(FIG, pattern = "\\.pdf$", full.names = TRUE)
file.copy(figs, file.path(OUT, "figures"), overwrite = TRUE)
cat("figures written:", basename(figs), "\n")
