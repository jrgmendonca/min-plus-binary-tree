#!/usr/bin/env Rscript
# Master figure regeneration script with improved aesthetics:
#   - Times-like (serif) font throughout
#   - Math-formatted log-axis tick labels (10^{-1}, 1, 10, 10^{2}, ...)
#     instead of "1e-01", "1e+00"
#   - Slightly thicker lines and bigger glyphs for visibility
#
# Usage:  Rscript redo_figures.R
# Outputs all six figures into figures/.

# ---- common style ----------------------------------------------------------

LWD  <- 1.6   # line width for data and reference lines
LWDR <- 1.2   # line width for reference (dashed) lines
PCH  <- 19
PCEX <- 0.95

set_style <- function() {
    par(family = "Times",
        font.lab = 1, font.axis = 1,
        cex.lab  = 1.20, cex.axis = 1.10,
        las  = 1,
        mgp  = c(2.8, 0.7, 0),
        mar  = c(4.5, 4.7, 0.8, 1.0),
        tcl  = -0.45)
}

log10_label <- function(p) {
    if (abs(p)     < 1e-9) quote(1)
    else if (abs(p - 1) < 1e-9) quote(10)
    else bquote(10^.(p))
}

draw_log10_axis <- function(side, lim_log10, minor = TRUE, ...) {
    major_p <- ceiling(lim_log10[1]):floor(lim_log10[2])
    major   <- 10^major_p
    labels  <- sapply(major_p, function(p) as.expression(log10_label(p)))
    axis(side, at = major, labels = labels, ...)
    if (minor) {
        minors <- as.vector(outer(2:9, major, "*"))
        keep   <- minors >= 10^lim_log10[1] & minors <= 10^lim_log10[2]
        axis(side, at = minors[keep], labels = FALSE, tcl = -0.25, ...)
    }
}

# ---- (a) critical_scaling --------------------------------------------------

{
    Ns <- c(20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60)
    files <- c(
        "runs/critical/N20_M100000_p0.5_ones.samples.dat",
        paste0("runs/critical/N", c(24, 28, 32), "_M100000_p0.5_K20.samples.dat"),
        paste0("runs/critical/N", c(36, 40, 44, 48, 52, 56, 60),
               "_M100000_p0.5_K35.samples.dat"))
    elog <- numeric(length(Ns)); selog <- numeric(length(Ns))
    for (i in seq_along(Ns)) {
        x <- scan(files[i], comment.char = "#", quiet = TRUE)
        lx <- log(x); elog[i] <- mean(lx); selog[i] <- sd(lx) / sqrt(length(lx))
    }
    slope_AC <- 2 * pi / (3 * sqrt(3))
    sqrtN    <- sqrt(Ns); inv_sqrtN <- 1 / sqrtN
    ratio    <- elog / sqrtN; ratio_se <- selog / sqrtN

    pdf("figures/critical_scaling.pdf", width = 7.5, height = 3.8, pointsize = 10)
    par(mfrow = c(1, 2),
        family = "Times",
        font.lab = 1, font.axis = 1,
        cex.lab = 1.20, cex.axis = 1.10,
        las = 1, mgp = c(2.8, 0.7, 0),
        mar = c(4.5, 4.7, 0.8, 0.8),
        tcl = -0.45)

    # Left
    plot(sqrtN, elog, type = "n",
         xlab = expression(sqrt(N)),
         ylab = expression(E*"["*log~X[N]*"]"),
         xlim = c(4, 8), ylim = c(0, 9))
    abline(0, slope_AC, lty = 2, col = "grey40", lwd = LWDR)
    arrows(sqrtN, elog - selog, sqrtN, elog + selog,
           length = 0.025, angle = 90, code = 3, lwd = LWD)
    points(sqrtN, elog, pch = PCH, cex = PCEX)
    text(6.5, 8.0, expression(slope == 2*pi/3*sqrt(3)),
         col = "grey40", cex = 0.95, pos = 2)
    grid(col = "grey85", lty = "dotted")

    # Right
    plot(inv_sqrtN, ratio, type = "n",
         xlab = expression(1/sqrt(N)),
         ylab = expression(E*"["*log~X[N]*"]"/sqrt(N)),
         xlim = c(0, 0.25), ylim = c(0.7, 1.25))
    abline(h = slope_AC, lty = 2, col = "grey40", lwd = LWDR)
    arrows(inv_sqrtN, ratio - ratio_se, inv_sqrtN, ratio + ratio_se,
           length = 0.025, angle = 90, code = 3, lwd = LWD)
    points(inv_sqrtN, ratio, pch = PCH, cex = PCEX)
    text(0, 1.235, expression(2*pi/3*sqrt(3)),
         col = "grey40", cex = 0.95, pos = 4)
    grid(col = "grey85", lty = "dotted")

    dev.off()
    cat("wrote figures/critical_scaling.pdf\n")
}

# ---- (b) qq_beta21 ---------------------------------------------------------

{
    N <- 60; cc <- pi^2 / 3
    s <- as.integer(scan("runs/critical/N60_M100000_p0.5_K35.samples.dat",
                         what = integer(), comment.char = "#", quiet = TRUE))
    ls <- log(s[s > 0])
    t  <- ls / sqrt(cc * N)
    probs <- seq(0.001, 0.999, length.out = 600)
    emp_q  <- quantile(t, probs = probs, type = 7)
    beta_q <- sqrt(probs)

    pdf("figures/qq_beta21.pdf", width = 4.5, height = 4.5, pointsize = 10)
    set_style()
    plot(beta_q, emp_q, type = "n",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = expression("Beta"~"(2,1) quantile"),
         ylab = expression(empirical~quantile~of~log(X[N])/sqrt(cN)))
    abline(0, 1, lty = 2, col = "grey40", lwd = LWDR)
    lines(beta_q, emp_q, col = "black", lwd = LWD)
    legend("topleft", inset = c(0.05, 0.04),
           c("empirical Q-Q at N = 60", "y = x  (AC limit)"),
           col = c("black", "grey40"), lty = c(1, 2),
           lwd = c(LWD, LWDR), bty = "n", cex = 1.10,
           y.intersp = 1.15)
    grid(col = "grey85", lty = "dotted")
    dev.off()
    cat("wrote figures/qq_beta21.pdf\n")
}

# ---- (c) dist_critical -----------------------------------------------------

{
    Ns_show <- c(20, 28, 36, 44, 52, 60)
    paths <- c(
        "runs/critical/N20_M100000_p0.5_ones.hist.csv",
        "runs/critical/N28_M100000_p0.5_K20.hist.csv",
        paste0("runs/critical/N", c(36, 44, 52, 60), "_M100000_p0.5_K35.hist.csv"))
    # Geometric (log-spaced) bin edges, deduplicated to integers
    edges <- unique(round(exp(seq(log(1), log(1500), length.out = 32))))
    bin_mid <- sqrt(head(edges, -1) * tail(edges, -1))
    bin_width <- diff(edges)
    # Two panels since the JPA revision: left, the distributions; right, the
    # mass-escape view at fixed k (Referee 2, Comment 2.3).
    pdf("figures/dist_critical.pdf", width = 7.5, height = 3.8, pointsize = 10)
    par(mfrow = c(1, 2),
        family = "Times",
        font.lab = 1, font.axis = 1,
        cex.lab  = 1.20, cex.axis = 1.10,
        las  = 1, mgp = c(2.8, 0.7, 0),
        mar  = c(4.5, 4.7, 0.8, 0.8),
        tcl  = -0.45)
    plot.new()
    plot.window(xlim = c(1, 1e3), ylim = c(1e-5, 2e-1), log = "xy")
    cols <- gray(seq(0.7, 0.0, length.out = length(Ns_show)))
    pchs <- c(15, 17, 19, 18, 4, 8)
    for (i in seq_along(Ns_show)) {
        d <- read.csv(paths[i], comment.char = "#")
        M_total <- sum(d$count)
        bin_idx <- findInterval(d$value, edges)
        valid <- bin_idx >= 1 & bin_idx <= length(edges) - 1 & d$count > 0
        if (!any(valid)) next
        bin_count <- tapply(d$count[valid], bin_idx[valid], sum)
        bins_used <- as.integer(names(bin_count))
        y <- as.numeric(bin_count) / M_total / bin_width[bins_used]
        x <- bin_mid[bins_used]
        keep <- y > 0 & x >= 1 & x <= 1e3
        points(x[keep], y[keep],
               col = cols[i], pch = pchs[i], cex = 0.95)
    }
    box()
    draw_log10_axis(1, c(0, 3))
    draw_log10_axis(2, c(-5, 0))
    title(xlab = expression(k), ylab = expression(widehat(P)*"("*X[N] == k*")"))
    legend("topright", inset = c(0.06, 0.02),
           legend = paste("N =", Ns_show),
           col = cols, pch = pchs, bty = "n", cex = 1.05,
           y.intersp = 1.15)
    grid(col = "grey85", lty = "dotted")

    # Right panel: P(X_N = k) vs 1/N at fixed k -- mass escapes to infinity.
    Ns_all <- c(20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60)
    paths_all <- c(
        "runs/critical/N20_M100000_p0.5_ones.hist.csv",
        paste0("runs/critical/N", c(24, 28, 32), "_M100000_p0.5_K20.hist.csv"),
        paste0("runs/critical/N", c(36, 40, 44, 48, 52, 56, 60),
               "_M100000_p0.5_K35.hist.csv"))
    kmax <- 5
    Pk <- matrix(0, nrow = length(Ns_all), ncol = kmax)
    for (i in seq_along(Ns_all)) {
        d <- read.csv(paths_all[i], comment.char = "#")
        Mtot <- sum(d$count)
        for (k in 1:kmax) {
            idx <- which(d$value == k)
            if (length(idx)) Pk[i, k] <- d$count[idx] / Mtot
        }
    }
    inv_N <- 1 / Ns_all
    colsk <- gray(seq(0.0, 0.62, length.out = kmax))
    pchsk <- c(19, 15, 17, 18, 4)
    plot(NULL, xlim = c(0, 0.054), ylim = c(0, 0.105),
         xlab = expression(1/N),
         ylab = expression(widehat(P)*"("*X[N] == k*")"))
    abline(0, 2, lty = 2, col = "grey40", lwd = LWDR)
    text(0.0465, 0.099, expression(2/N), col = "grey40", cex = 0.95, pos = 2)
    for (k in 1:kmax) {
        lines(inv_N, Pk[, k], lwd = 1.2, col = colsk[k])
        points(inv_N, Pk[, k], pch = pchsk[k], cex = 0.85, col = colsk[k])
    }
    legend("topleft", inset = c(0.02, 0.02),
           legend = c(expression(k == 1), expression(k == 2),
                      expression(k == 3), expression(k == 4),
                      expression(k == 5)),
           col = colsk, pch = pchsk, lty = 1, lwd = 1.2,
           bty = "n", cex = 1.0, y.intersp = 1.1)
    grid(col = "grey85", lty = "dotted")

    dev.off()
    cat("wrote figures/dist_critical.pdf\n")
}

# ---- (d) noncritical -------------------------------------------------------

{
    ps <- c(0.10, 0.20, 0.30, 0.40)
    paths <- sprintf("runs/noncritical/N48_M100000_p%.2f_K35.hist.csv", ps)

    pdf("figures/noncritical.pdf", width = 5.5, height = 4.5, pointsize = 10)
    set_style()
    par(cex.lab = 1.30, cex.axis = 1.20,
        mgp = c(2.9, 0.75, 0), mar = c(4.7, 4.9, 0.8, 1.0))
    plot.new()
    plot.window(xlim = c(1, 12), ylim = c(1e-4, 1e0), log = "y")
    cols <- gray(seq(0.0, 0.55, length.out = length(ps)))
    pchs <- c(15, 17, 19, 18)
    ltys <- c(2, 3, 4, 5)
    for (i in seq_along(ps)) {
        d <- read.csv(paths[i], comment.char = "#")
        d$prob <- d$count / sum(d$count)
        keep <- d$value > 0 & d$prob > 0 & d$value <= 12
        lines(d$value[keep], d$prob[keep],
              col = cols[i], lwd = LWD, lty = ltys[i])
        points(d$value[keep], d$prob[keep],
               col = cols[i], pch = pchs[i], cex = PCEX)
    }
    axis(1, at = seq(2, 12, by = 2))
    draw_log10_axis(2, c(-4, 0))
    box()
    title(xlab = expression(k), ylab = expression(widehat(P)*"("*X[N] == k*")"))
    legend("topright", inset = c(0.05, 0.04),
           legend = c(expression(p == 0.1), expression(p == 0.2),
                      expression(p == 0.3), expression(p == 0.4)),
           col = cols, pch = pchs, lty = ltys,
           lwd = LWD, bty = "n", cex = 1.10,
           y.intersp = 1.15)
    grid(col = "grey85", lty = "dotted")
    dev.off()
    cat("wrote figures/noncritical.pdf\n")
}

# ---- (e) bern_p0 -----------------------------------------------------------

{
    qs <- sprintf("%.2f", seq(0, 1, by = 0.05))
    p0 <- numeric(length(qs))
    for (i in seq_along(qs)) {
        f <- sprintf("runs/bernoulli/N48_M100000_q%s_K35.summary.txt", qs[i])
        ln <- readLines(f)
        k0 <- ln[grepl("^\\s*k=0\\s+", ln)][1]
        p0[i] <- as.numeric(sub(".*prob=", "", k0))
    }
    qv <- as.numeric(qs)

    pdf("figures/bern_p0.pdf", width = 4.8, height = 4.5, pointsize = 10)
    set_style()
    plot(qv, p0, type = "n",
         xlim = c(0, 1), ylim = c(0, 1),
         xlab = expression(q),
         ylab = expression(widehat(P)*"("*X[N] == 0*")"))
    abline(0, 1, lty = 2, col = "grey40", lwd = LWDR)
    points(qv, p0, pch = PCH, cex = PCEX + 0.1)
    legend("topleft", inset = c(0.05, 0.04),
           c("data", expression(P*"("*X[N] == 0*") =" ~ q)),
           col = c("black", "grey40"), pch = c(PCH, NA),
           lty = c(NA, 2), lwd = c(NA, LWDR),
           bty = "n", cex = 1.10,
           y.intersp = 1.15)
    grid(col = "grey85", lty = "dotted")
    dev.off()
    cat("wrote figures/bern_p0.pdf\n")
}

# ---- (f) bern_moments ------------------------------------------------------

{
    qs <- sprintf("%.2f", seq(0, 1, by = 0.05))
    qv <- as.numeric(qs)
    means <- vars <- numeric(length(qs))
    for (i in seq_along(qs)) {
        f <- sprintf("runs/bernoulli/N48_M100000_q%s_K35.summary.txt", qs[i])
        ln <- readLines(f)
        means[i] <- as.numeric(sub(".*:", "", ln[grepl("mean\\(X_N\\)", ln)][1]))
        vars[i]  <- as.numeric(sub(".*:", "", ln[grepl("var\\(X_N\\)",  ln)][1]))
    }
    keep <- means > 0 & vars > 0   # drop q = 1 (both zero)

    pdf("figures/bern_moments.pdf", width = 5.4, height = 4.5, pointsize = 10)
    par(family = "Times",
        font.lab = 1, font.axis = 1,
        cex.lab = 1.30, cex.axis = 1.20,
        las = 1, mgp = c(2.9, 0.75, 0),
        mar = c(4.7, 4.9, 0.8, 5.0),
        tcl = -0.45)
    # Mean (left axis)
    plot(qv[keep], means[keep], type = "n", log = "y", yaxt = "n",
         xlim = c(0, 1), ylim = c(1e-1, 1e4),
         xlab = expression(q),
         ylab = expression(widehat(E)*"["*X[N]*"]"))
    draw_log10_axis(2, c(-1, 4))
    points(qv[keep], means[keep], pch = PCH, cex = PCEX)
    lines(qv[keep], means[keep], lwd = LWD)
    # Variance (right axis)
    par(new = TRUE)
    plot(qv[keep], vars[keep], type = "n", log = "y",
         xlim = c(0, 1), ylim = c(1e0, 1e7),
         axes = FALSE, xlab = "", ylab = "")
    draw_log10_axis(4, c(0, 7), col.axis = "grey40")
    mtext(expression(widehat(Var)(X[N])), side = 4, line = 2.8,
          las = 0, col = "grey40", cex = 1.05)
    points(qv[keep], vars[keep], pch = 1, cex = PCEX, col = "grey40")
    lines(qv[keep], vars[keep], lwd = LWD, col = "grey40", lty = 3)

    legend("topright", inset = c(0.05, 0.04),
           c("mean", "variance"),
           col = c("black", "grey40"), pch = c(PCH, 1),
           lty = c(1, 3), lwd = LWD, bty = "n", cex = 1.10,
           y.intersp = 1.15)
    grid(col = "grey85", lty = "dotted")
    dev.off()
    cat("wrote figures/bern_moments.pdf\n")
}

cat("\nAll figures regenerated.\n")
