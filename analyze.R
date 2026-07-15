# analyze.R - load and analyse min-plus simulator output.
#
# Interactive use:
#   source("analyze.R")
#   runs <- discover_runs("runs")
#   s <- summary_table(runs); print(s)
#   plot_bern_transition(runs, depth = 26)
#   plot_critical_scaling(runs)
#   plot_distribution("runs/critical/N26_M10000_p0.5_four")
#
# Batch use (Rscript):
#   Rscript analyze.R [outdir]      # writes summary.csv + canonical pdfs
#
# Reads only the histogram CSVs by default; samples.dat is loaded on demand
# via load_samples(prefix). Base R only -- no external packages.

# ---------- I/O ----------

# Parse the metadata header of a hist.csv file. Header line 2 looks like:
#   # N=22  p=0.500..  q=0.300..  M=5000  leaves=3  seed=0x...  wall=6.423s
parse_hist_header <- function(path) {
    hdr <- readLines(path, n = 2)
    if (length(hdr) < 2) stop("missing header in ", path)
    fields <- strsplit(sub("^#\\s*", "", hdr[2]), "\\s+")[[1]]
    out <- list(path = path)
    for (f in fields) {
        kv <- strsplit(f, "=", fixed = TRUE)[[1]]
        if (length(kv) == 2) {
            key <- kv[1]
            raw <- kv[2]
            num <- suppressWarnings(as.numeric(sub("s$", "", raw)))
            out[[key]] <- if (!is.na(num)) num else raw
        }
    }
    out
}

# Read a histogram CSV into a data.frame {value, count}; metadata in attr.
load_hist <- function(prefix) {
    path <- paste0(prefix, ".hist.csv")
    h <- read.csv(path, comment.char = "#")
    attr(h, "meta") <- parse_hist_header(path)
    h
}

# Read raw samples (one integer per line). Larger than load_hist; lazy.
load_samples <- function(prefix) {
    scan(paste0(prefix, ".samples.dat"),
         what = integer(), comment.char = "#", quiet = TRUE)
}

# Recursively find all run prefixes under a directory.
discover_runs <- function(root = "runs") {
    files <- list.files(root, pattern = "\\.hist\\.csv$",
                        recursive = TRUE, full.names = TRUE)
    sub("\\.hist\\.csv$", "", files)
}

# ---------- summary ----------

# Stats from a histogram (avoids re-reading samples.dat).
hist_stats <- function(h) {
    n <- sum(h$count)
    mean_x <- sum(h$value * h$count) / n
    var_x  <- sum(h$count * (h$value - mean_x)^2) / n
    p_at <- function(k) {
        idx <- which(h$value == k)
        if (length(idx)) h$count[idx] / n else 0
    }
    list(M = n, mean = mean_x, var = var_x, sd = sqrt(var_x),
         min = min(h$value[h$count > 0]),
         max = max(h$value[h$count > 0]),
         p0 = p_at(0), p1 = p_at(1), p2 = p_at(2))
}

# Build a wide summary across runs.
summary_table <- function(prefixes) {
    rows <- lapply(prefixes, function(p) {
        meta <- parse_hist_header(paste0(p, ".hist.csv"))
        h <- read.csv(paste0(p, ".hist.csv"), comment.char = "#")
        st <- hist_stats(h)
        cls <- basename(dirname(p))
        leaves_code <- meta$leaves %||% NA
        # effective depth: leaves=2 (four) shifts +2; leaves=4 (bern-shift1) +1
        eff <- meta$N + switch(as.character(leaves_code),
                               "0" = 0, "1" = 1, "2" = 2,
                               "3" = 0, "4" = 1, NA)
        data.frame(prefix = basename(p), class = cls,
                   N = meta$N, eff_N = eff, p = meta$p, q = meta$q %||% NA,
                   leaves = leaves_code, M = st$M,
                   mean = st$mean, var = st$var, sd = st$sd,
                   p0 = st$p0, p1 = st$p1, p2 = st$p2,
                   max = st$max,
                   stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------- plots (base R only) ----------

plot_bern_transition <- function(prefixes, depth = NULL, outpdf = NULL) {
    s <- summary_table(prefixes[grepl("/bernoulli/", prefixes)])
    if (nrow(s) == 0) { message("no bernoulli runs"); return(invisible(NULL)) }
    if (!is.null(depth)) s <- s[s$N == depth, ]
    s <- s[order(s$q), ]
    if (!is.null(outpdf)) pdf(outpdf, width = 5.5, height = 4)
    plot(s$q, s$p0, type = "b", pch = 19, ylim = c(0, 1),
         xlab = expression(q), ylab = expression(P(X[N] == 0)),
         main = sprintf("Bernoulli variant, N=%d, M=%d", s$N[1], s$M[1]))
    abline(h = 0.5, lty = 3, col = "grey50")
    grid()
    if (!is.null(outpdf)) { dev.off(); message("wrote ", outpdf) }
    invisible(s)
}

plot_critical_scaling <- function(prefixes, outpdf = NULL) {
    s <- summary_table(prefixes[grepl("/critical/", prefixes)])
    if (nrow(s) == 0) { message("no critical runs"); return(invisible(NULL)) }
    s <- s[order(s$eff_N), ]
    if (!is.null(outpdf)) pdf(outpdf, width = 8, height = 4)
    op <- par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    on.exit(par(op))
    plot(s$eff_N, s$mean, type = "b", pch = 19, log = "y",
         xlab = "effective depth", ylab = expression(E*"["*X[N]*"]"),
         main = "mean")
    grid()
    plot(s$eff_N, s$var, type = "b", pch = 19, log = "y",
         xlab = "effective depth", ylab = expression(Var(X[N])),
         main = "variance")
    grid()
    if (!is.null(outpdf)) { dev.off(); message("wrote ", outpdf) }
    invisible(s)
}

plot_distribution <- function(prefix, outpdf = NULL,
                              log = "y", xmax = NULL) {
    h <- load_hist(prefix)
    meta <- attr(h, "meta")
    n <- sum(h$count); h$prob <- h$count / n
    if (!is.null(xmax)) h <- h[h$value <= xmax, ]
    if (!is.null(outpdf)) pdf(outpdf, width = 5.5, height = 4)
    title <- sprintf("X_N distribution (N=%g, p=%g, q=%g, M=%g)",
                     meta$N, meta$p, meta$q %||% 0, meta$M)
    plot(h$value, h$prob, type = "h", xlab = expression(X[N]),
         ylab = "probability", log = log, main = title)
    if (!is.null(outpdf)) { dev.off(); message("wrote ", outpdf) }
    invisible(h)
}

# ---------- batch entrypoint ----------

if (sys.nframe() == 0L) {
    args <- commandArgs(trailingOnly = TRUE)
    outdir <- if (length(args) >= 1) args[1] else "figures"
    dir.create(outdir, showWarnings = FALSE)

    runs <- discover_runs("runs")
    if (length(runs) == 0L) { message("no runs under runs/"); quit(status = 1) }

    s <- summary_table(runs)
    write.csv(s, file.path(outdir, "summary.csv"), row.names = FALSE)
    message("wrote ", file.path(outdir, "summary.csv"))
    print(s, row.names = FALSE)

    if (any(s$class == "bernoulli")) {
        for (n in unique(s$N[s$class == "bernoulli"])) {
            plot_bern_transition(runs, depth = n,
                outpdf = file.path(outdir, sprintf("bern_transition_N%d.pdf", n)))
        }
    }
    if (any(s$class == "critical")) {
        plot_critical_scaling(runs,
            outpdf = file.path(outdir, "critical_scaling.pdf"))
        # one distribution sample at the largest available eff_N
        sc <- s[s$class == "critical", ]
        biggest <- sc$prefix[which.max(sc$eff_N)]
        path <- file.path("runs/critical", biggest)
        plot_distribution(path,
            outpdf = file.path(outdir, sprintf("dist_%s.pdf", biggest)))
    }
}
