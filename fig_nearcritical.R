#!/usr/bin/env Rscript
# fig_nearcritical.R - two-panel figure for the near-critical scaling window
# p_N = 1/2 - c/N (JPA revision, new subsection 3.3).
#   Left : A(c,N) = E[log X_N]/sqrt(N) vs 1/sqrt(N) for c = 0, 1/2, 1, 2,
#          with the critical asymptote 2*pi/(3*sqrt(3)) dashed.
#   Right: deficit phi(c,N) = A(0,N) - A(c,N), nearly constant in N.
# Style matches redo_figures.R (Times family, grayscale, dotted grid).
# Usage: Rscript fig_nearcritical.R   (from the simulator directory)

LWD  <- 1.6
LWDR <- 1.2
PCEX <- 0.95

elog_from_hist <- function(path) {
    d <- read.csv(path, comment.char = "#")
    M <- sum(d$count)
    lv <- log(d$value)
    m  <- sum(d$count * lv) / M
    v  <- sum(d$count * lv^2) / M - m^2
    c(mean = m, se = sqrt(v / M))
}

Ns <- c(24, 32, 40, 48, 56, 60)
Ks <- ifelse(Ns <= 32, 20, 35)

paths <- list(
    "0"   = sprintf("runs/critical/N%d_M100000_p0.5_K%d.hist.csv", Ns, Ks),
    "0.5" = c("runs/nearcritical/N24_M100000_p0.4791666667_c05_K20.hist.csv",
              "runs/nearcritical/N32_M100000_p0.484375_c05_K20.hist.csv",
              "runs/nearcritical/N40_M100000_p0.4875_c05_K35.hist.csv",
              "runs/nearcritical/N48_M100000_p0.4895833333_c05_K35.hist.csv",
              "runs/nearcritical/N56_M100000_p0.4910714286_c05_K35.hist.csv",
              "runs/nearcritical/N60_M100000_p0.4916666667_c05_K35.hist.csv"),
    "1"   = c("runs/nearcritical/N24_M100000_p0.4583333333_c10_K20.hist.csv",
              "runs/nearcritical/N32_M100000_p0.46875_c10_K20.hist.csv",
              "runs/nearcritical/N40_M100000_p0.475_c10_K35.hist.csv",
              "runs/nearcritical/N48_M100000_p0.4791666667_c10_K35.hist.csv",
              "runs/nearcritical/N56_M100000_p0.4821428571_c10_K35.hist.csv",
              "runs/nearcritical/N60_M100000_p0.4833333333_c10_K35.hist.csv"),
    "2"   = c("runs/nearcritical/N24_M100000_p0.4166666667_c20_K20.hist.csv",
              "runs/nearcritical/N32_M100000_p0.4375_c20_K20.hist.csv",
              "runs/nearcritical/N40_M100000_p0.45_c20_K35.hist.csv",
              "runs/nearcritical/N48_M100000_p0.4583333333_c20_K35.hist.csv",
              "runs/nearcritical/N56_M100000_p0.4642857143_c20_K35.hist.csv",
              "runs/nearcritical/N60_M100000_p0.4666666667_c20_K35.hist.csv"))

A  <- list()
SE <- list()
for (cc in names(paths)) {
    st <- sapply(paths[[cc]], elog_from_hist)
    A[[cc]]  <- st["mean", ] / sqrt(Ns)
    SE[[cc]] <- st["se", ]   / sqrt(Ns)
}

x <- 1 / sqrt(Ns)
slope_AC <- 2 * pi / (3 * sqrt(3))
cols <- c(gray(0), gray(0.30), gray(0.50), gray(0.62))
pchs <- c(19, 15, 17, 18)
ltys <- c(1, 2, 4, 5)
cexs <- c(PCEX, PCEX, PCEX, PCEX + 0.15)   # pch 18 renders smaller
cset <- names(paths)

pdf("figures/nearcritical.pdf", width = 7.5, height = 3.8, pointsize = 10)
par(mfrow = c(1, 2),
    family = "Times",
    font.lab = 1, font.axis = 1,
    cex.lab = 1.20, cex.axis = 1.10,
    las = 1, mgp = c(2.8, 0.7, 0),
    mar = c(4.5, 4.7, 0.8, 0.8),
    tcl = -0.45)

# Left: A(c, N) vs 1/sqrt(N)
plot(NULL, xlim = c(0, 0.22), ylim = c(0.25, 1.25),
     xlab = expression(1/sqrt(N)),
     ylab = expression(widehat(A)*"("*c*","~N*")" == widehat(E)*"["*log~X[N]*"]"/sqrt(N)))
abline(h = slope_AC, lty = 2, col = "grey40", lwd = LWDR)
text(0, 1.243, expression(2*pi/3*sqrt(3)), col = "grey40", cex = 0.95, pos = 4)
for (i in seq_along(cset)) {
    a <- A[[cset[i]]]; s <- SE[[cset[i]]]
    lines(x, a, lty = ltys[i], lwd = 1.3, col = cols[i])
    arrows(x, a - s, x, a + s, length = 0.02, angle = 90, code = 3,
           lwd = LWD, col = cols[i])
    points(x, a, pch = pchs[i], cex = cexs[i], col = cols[i])
}
legend("bottomleft", inset = c(0.02, 0.03),
       legend = c(expression(c == 0 ~ "(critical)"), expression(c == 1/2),
                  expression(c == 1), expression(c == 2)),
       col = cols, pch = pchs, lty = ltys, lwd = 1.3,
       bty = "n", cex = 1.02, y.intersp = 1.15)
grid(col = "grey85", lty = "dotted")

# Right: deficit phi(c, N) = A(0, N) - A(c, N)
plot(NULL, xlim = c(0, 0.22), ylim = c(0.10, 0.55),
     xlab = expression(1/sqrt(N)),
     ylab = expression(widehat(varphi)*"("*c*","~N*")" ==
                       widehat(A)*"("*0*","~N*")" - widehat(A)*"("*c*","~N*")"))
for (i in 2:4) {
    ph  <- A[["0"]] - A[[cset[i]]]
    phs <- sqrt(SE[["0"]]^2 + SE[[cset[i]]]^2)
    lines(x, ph, lty = ltys[i], lwd = 1.3, col = cols[i])
    arrows(x, ph - phs, x, ph + phs, length = 0.02, angle = 90, code = 3,
           lwd = LWD, col = cols[i])
    points(x, ph, pch = pchs[i], cex = cexs[i], col = cols[i])
    text(0.002, mean(ph), bquote(c == .(c("", "1/2", "1", "2")[i])),
         pos = 4, cex = 1.02)
}
grid(col = "grey85", lty = "dotted")

dev.off()
cat("wrote figures/nearcritical.pdf\n")

# Console table for the manuscript
tab <- do.call(cbind, A)
rownames(tab) <- Ns
cat("\nA(c, N):\n"); print(round(tab, 4))
cat("\nphi(c, N):\n"); print(round(tab[, 1] - tab[, 2:4], 4))
