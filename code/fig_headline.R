###############################################################################
## Headline figure (Figure 1) of
##   "What does a ridge network estimate? Effective penalty and leakage in
##    ridge precision matrix estimation"
##
## Author : M. Arashi
## Purpose: one two-panel figure for the start of Section 3, both panels for
##          the band model of Table 2 (p = 200, huge band generator, which
##          uses no random numbers).
##          (a) effective penalty against nominal penalty: the root kappa-bar
##              of the fixed-point equation (9) for n = 50 and n = 400, the
##              identity, the Ky Fan bound of (8) for n = 50 and the Monte
##              Carlo kappa_n of Table 2;
##          (b) leakage: mean |rho^(kappa)_ij| of the population ridge network
##              by graph distance, with the bound of Proposition 3.4(c).
##          Everything except the Monte Carlo points is deterministic.
## Input  : ../results/sim_kappa_grid.csv (Monte Carlo kappa_n)
## Output : ../figs/fig_headline.pdf, copied to ../figures/
## Usage  : Rscript fig_headline.R   (run from code/); also called at the end
##          of make_figures.R
###############################################################################

lib_extra <- Sys.getenv("RIDGE_RLIB")
if (nzchar(lib_extra)) .libPaths(c(lib_extra, .libPaths()))
suppressPackageStartupMessages(library(huge))
OUT <- ".."
FIG <- file.path(OUT, "figs")
dir.create(FIG, showWarnings = FALSE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE)

okabe <- c(orange = "#E69F00", skyblue = "#56B4E9", green = "#009E73",
           blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7",
           black = "#000000")

## band model of Table 2
p <- 200
Sigma <- as.matrix(huge.generator(d = p, graph = "band", g = 1, v = 0.4, u = 0.1,
                                  verbose = FALSE)$sigma)
ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
smax <- max(ev); tmax <- 1 / min(ev); tmin <- 1 / smax

## root of (9) on (lambda, infinity)
kappa_bar <- function(n, lambda) {
  psi <- function(k) lambda / k + sum(ev / (ev + k)) / n - 1
  uniroot(psi, c(lambda, lambda + sum(ev) / n + 1e-12), tol = 1e-14)$root
}
## lambda-free part of the Ky Fan bound of (8)
kyfan <- function(n) {
  k <- max(0, p - n + 1)
  if (k > 0) sum(sort(ev)[1:k]) / n else 0
}

lam <- 10^seq(-3, 1, length.out = 161)
kb50 <- sapply(lam, kappa_bar, n = 50)
kb400 <- sapply(lam, kappa_bar, n = 400)
kf50 <- kyfan(50)
kg <- read.csv(file.path(OUT, "results", "sim_kappa_grid.csv"))

## population ridge partial correlations by graph distance (a chain, so the
## distance of i and j is |i - j|)
kaps <- c(0.1, 1, 3, 10)
dmax <- 6
D <- abs(outer(1:p, 1:p, "-"))
pc_true <- -cov2cor(solve(Sigma))[1, 2]
leak <- sapply(kaps, function(k) {
  R <- -cov2cor(solve(Sigma + k * diag(p)))
  sapply(1:dmax, function(d) mean(abs(R[D == d])))
})
bound <- sapply(kaps, function(k) {
  vr <- (1 + k * tmax) / (1 + k * tmin)
  q <- (sqrt(vr) - 1) / (sqrt(vr) + 1)
  smax * (1 + sqrt(vr))^2 / (2 * k * vr) * q^(1:dmax)
})

## labels 10^k on a log axis
pow10_axis <- function(side, k) {
  axis(side, at = 10^k, labels = as.expression(lapply(k, function(e) bquote(10^.(e)))),
       las = 1)
}

## drawn at the printed size (text width 5.6 in), so 9 pt here is 9 pt on paper
cairo_pdf(file.path(FIG, "fig_headline.pdf"), width = 5.6, height = 2.8,
          pointsize = 9, family = "sans")
par(mfrow = c(1, 2), mar = c(3.3, 3.9, 1.3, 0.6), oma = c(0, 0, 0, 0),
    mgp = c(2.2, 0.55, 0), tcl = -0.35, las = 1, cex = 1)

## (a) effective penalty
c50 <- okabe["vermillion"]; c400 <- okabe["blue"]
plot(NA, xlim = c(1e-3, 10), ylim = c(1e-3, 20), log = "xy", xaxt = "n", yaxt = "n",
     xlab = expression("nominal penalty " * lambda),
     ylab = expression("effective penalty " * kappa))
pow10_axis(1, -3:1); pow10_axis(2, -3:1)
lines(lam, lam, col = "gray55", lwd = 1.2)
lines(lam, lam + kf50, col = c50, lty = 2, lwd = 1.3)
lines(lam, kb400, col = c400, lwd = 2)
lines(lam, kb50, col = c50, lwd = 2)
m50 <- kg[kg$n == 50, ]; m400 <- kg[kg$n == 400, ]
points(m50$lambda, m50$kappa_mc, pch = 21, bg = "white", col = c50, cex = 1.1, lwd = 1.3)
points(m400$lambda, m400$kappa_mc, pch = 22, bg = "white", col = c400, cex = 1.05, lwd = 1.3)
text(1.2e-3, kb50[1] * 1.45, expression(paste(italic(n) == 50, "  (", italic(p) / italic(n) == 4, ")")),
     col = c50, adj = c(0, 0))
text(1.2e-3, kf50 / 1.45, "Ky Fan bound, n = 50", col = c50, adj = c(0, 1), font = 1)
text(2.5e-2, 2.5e-2 / 3.2,
     expression(paste(italic(n) == 400, "  (", italic(p) / italic(n) == 1/2, ")")),
     col = c400, adj = c(0, 1))
text(2.2, 2.2 / 2.2, expression(kappa == lambda), col = "gray40", adj = c(0, 1))
## key for the Monte Carlo points: both symbols, then the words
points(c(7e-3, 1.12e-2), c(1.9e-3, 1.9e-3), pch = c(21, 22), col = c(c50, c400),
       bg = "white", cex = 1.05, lwd = 1.3)
text(1.6e-2, 1.9e-3, expression("Monte Carlo " * kappa[italic(n)]), adj = c(0, 0.5))
mtext("(a)", side = 3, line = 0.25, adj = 0, cex = 1, font = 2)

## (b) leakage by graph distance
kcol <- okabe[c("skyblue", "blue", "vermillion", "orange")]
kpch <- c(16, 15, 17, 18)
plot(NA, xlim = c(0.8, dmax + 0.2), ylim = c(1e-7, 1), log = "y", xaxt = "n", yaxt = "n",
     xlab = expression("graph distance " * italic(d)),
     ylab = expression("mean |" * rho[italic(ij)]^{(kappa)} * "|"))
axis(1, at = 1:dmax)
pow10_axis(2, seq(-7, 0, by = 1))
for (j in seq_along(kaps)) {
  lines(2:dmax, bound[2:dmax, j], col = kcol[j], lty = 3, lwd = 1.1)
  lines(1:dmax, leak[, j], col = kcol[j], lwd = 1.6)
  points(1:dmax, leak[, j], col = kcol[j], pch = kpch[j], cex = if (kpch[j] == 18) 1.35 else 1)
}
## the partial correlations of Theta: 0.40 on edges (d = 1), zero beyond
points(0.86, abs(pc_true), pch = 8, cex = 1.2, lwd = 1.3)
text(1.12, abs(pc_true) * 1.15, bquote(Theta * ": " * .(sprintf("%.2f", abs(pc_true))) * " on edges, 0 beyond"), adj = c(-0.08, 0.5))
legend("bottomleft", bty = "n", cex = 1, y.intersp = 1,
       legend = as.expression(c(lapply(kaps, function(k) bquote(kappa == .(k))),
                                "bound, Prop. 3.4(c)")),
       col = c(kcol, "gray30"), pch = c(kpch, NA), lty = c(1, 1, 1, 1, 3),
       lwd = c(1.6, 1.6, 1.6, 1.6, 1.1), pt.cex = c(1, 1, 1, 1.35, 1))
mtext("(b)", side = 3, line = 0.25, adj = 0, cex = 1, font = 2)
dev.off()

file.copy(file.path(FIG, "fig_headline.pdf"), file.path(OUT, "figures"), overwrite = TRUE)

## numbers quoted in the text and caption
cat(sprintf("Ky Fan bound (lambda-free), n = 50: %.3f\n", kf50))
cat(sprintf("kappa-bar at lambda = 1e-3: n = 50 %.3f, n = 400 %.5f\n", kb50[1], kb400[1]))
cat(sprintf("kappa-bar at lambda = 10:   n = 50 %.3f, n = 400 %.3f\n",
            kb50[length(lam)], kb400[length(lam)]))
cat(sprintf("true |partial correlation| of an edge: %.3f\n", abs(pc_true)))
print(data.frame(kappa = kaps, edge = leak[1, ], d2 = leak[2, ], d6 = leak[dmax, ],
                 rate = leak[3, ] / leak[2, ]), digits = 3)
