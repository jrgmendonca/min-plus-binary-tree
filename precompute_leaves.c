/*
 * precompute_leaves.c -- exact level-k distribution of the min-plus
 * recursion, with FFT-based autoconvolution for large k.
 *
 * Iterates
 *   mu_{k+1}(j) = p * (mu_k * mu_k)(j)
 *               + (1-p) * [ barF_k(j)^2 - barF_k(j+1)^2 ]
 * starting from the Bernoulli initial condition mu_0(0) = q, mu_0(1) = 1-q,
 * where (mu * mu)(j) = sum_{i} mu(i) mu(j-i) is the autoconvolution and
 * barF_k(j) = P(X_k >= j) is the right-tail cdf. At k = 1, 2 this reproduces
 * the one- and two-level leaf constructions of the companion paper (Appendix
 * A.1) and it matches the y-recursion of Lemma 1 along the j = 1 slice.
 *
 * The basic algorithm is O(4^K / 3) doubles of work and is fine up to
 * K ~ 18 (~30 s). For larger K the support is capped at 2^M (default
 * M = 20), discarding deep-tail mass beyond j_max = 2^M, and the
 * autoconvolution is computed via a radix-2 Cooley-Tukey FFT in
 * O(K * 2^M log 2^M) doubles. The truncated mass per level is reported
 * to stderr; for the critical case at K = 35 with M = 20 it is < 10^-15.
 *
 * Output is a CSV with one (value, prob) row per nonzero entry plus a
 * metadata header.
 *
 * Build:  make precompute_leaves
 * Usage:  precompute_leaves --prob P [--q Q] --levels K --out FILE
 *                           [--buffer-bits M] [--no-fft]
 */
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- radix-2 Cooley-Tukey FFT ---------- */

typedef struct { double re, im; } cplx;

/* In-place radix-2 FFT. n must be a power of 2.
 * sign = -1 forward, +1 inverse (caller normalises). */
static void fft_inplace(cplx *a, size_t n, int sign) {
    /* Bit-reversal permutation */
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        while (j & bit) { j ^= bit; bit >>= 1; }
        j ^= bit;
        if (i < j) { cplx t = a[i]; a[i] = a[j]; a[j] = t; }
    }
    /* Butterflies */
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = sign * 2.0 * M_PI / (double)len;
        const double wnre = cos(ang), wnim = sin(ang);
        const size_t half = len >> 1;
        for (size_t i = 0; i < n; i += len) {
            double wre = 1.0, wim = 0.0;
            for (size_t k = 0; k < half; k++) {
                const cplx u = a[i + k];
                const double bre = a[i + k + half].re;
                const double bim = a[i + k + half].im;
                const double vre = bre * wre - bim * wim;
                const double vim = bre * wim + bim * wre;
                a[i + k].re = u.re + vre;
                a[i + k].im = u.im + vim;
                a[i + k + half].re = u.re - vre;
                a[i + k + half].im = u.im - vim;
                const double nwre = wre * wnre - wim * wnim;
                const double nwim = wre * wnim + wim * wnre;
                wre = nwre; wim = nwim;
            }
        }
    }
}

/* Real autoconvolution via complex FFT.
 *   in[0..n_in-1] is real;
 *   out[0..n_out-1] receives (in * in)[0..n_out-1];
 *   scratch is a complex buffer of length n_fft (must be a power of 2,
 *   n_fft >= 2*n_in to avoid circular wraparound).
 */
static void real_autoconv_fft(const double * __restrict__ in, size_t n_in,
                              double * __restrict__ out, size_t n_out,
                              cplx * __restrict__ scratch, size_t n_fft) {
    for (size_t i = 0; i < n_in; i++) {
        scratch[i].re = in[i]; scratch[i].im = 0.0;
    }
    for (size_t i = n_in; i < n_fft; i++) {
        scratch[i].re = 0.0;  scratch[i].im = 0.0;
    }
    fft_inplace(scratch, n_fft, -1);
    for (size_t i = 0; i < n_fft; i++) {
        const double r = scratch[i].re, im = scratch[i].im;
        scratch[i].re = r*r - im*im;
        scratch[i].im = 2.0 * r * im;
    }
    fft_inplace(scratch, n_fft, +1);
    const double inv_n = 1.0 / (double)n_fft;
    for (size_t i = 0; i < n_out; i++) {
        double v = scratch[i].re * inv_n;
        if (v < 0.0) v = 0.0;   /* clamp roundoff noise */
        out[i] = v;
    }
}

/* Basic O(N^2) autoconvolution, used at small support. */
static void real_autoconv_basic(const double * __restrict__ in, size_t n_in,
                                double * __restrict__ out, size_t n_out) {
    for (size_t j = 0; j < n_out; j++) {
        double s = 0.0;
        const size_t lo = (j + 1 > n_in) ? j + 1 - n_in : 0;
        const size_t hi = (j < n_in) ? j : n_in - 1;
        for (size_t i = lo; i <= hi; i++) s += in[i] * in[j - i];
        out[j] = s;
    }
}

/* ---------- main recursion ---------- */

static void usage(const char *p) {
    fprintf(stderr,
        "usage: %s --prob P [--q Q] --levels K --out FILE\n"
        "          [--buffer-bits M] [--no-fft]\n"
        "  Compute mu_k(j) = P(X_k = j) by iterating the min-plus recursion\n"
        "  from the Bernoulli initial condition mu_0(0)=q, mu_0(1)=1-q.\n"
        "  Buffer caps the support at 2^M values (default M=20). FFT-based\n"
        "  autoconvolution kicks in once support exceeds 2^12.\n", p);
    exit(2);
}

int main(int argc, char **argv) {
    double p = -1.0, q = 0.0;
    int K = -1, buffer_bits = 20, no_fft = 0;
    const char *out_path = NULL;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if      (!strcmp(a, "--prob") && i+1 < argc)        p = atof(argv[++i]);
        else if (!strcmp(a, "--q") && i+1 < argc)           q = atof(argv[++i]);
        else if (!strcmp(a, "--levels") && i+1 < argc)      K = atoi(argv[++i]);
        else if (!strcmp(a, "--buffer-bits") && i+1 < argc) buffer_bits = atoi(argv[++i]);
        else if (!strcmp(a, "--no-fft"))                    no_fft = 1;
        else if (!strcmp(a, "--out") && i+1 < argc)         out_path = argv[++i];
        else usage(argv[0]);
    }
    if (p < 0.0 || p > 1.0 || q < 0.0 || q > 1.0 || K < 0 || K > 60 ||
        buffer_bits < 4 || buffer_bits > 25 || !out_path)
        usage(argv[0]);

    const size_t cap   = (size_t)1 << buffer_bits;          /* support <= cap */
    const size_t fft_n = (size_t)1 << (buffer_bits + 1);    /* 2 * cap */

    double *mu       = (double *)calloc(cap,         sizeof(double));
    double *mu_new   = (double *)calloc(cap,         sizeof(double));
    double *barF     = (double *)calloc(cap + 1,     sizeof(double));
    double *conv_buf = (double *)calloc(cap,         sizeof(double));
    cplx   *fft_scr  = (cplx   *)calloc(fft_n,       sizeof(cplx));
    if (!mu || !mu_new || !barF || !conv_buf || !fft_scr) {
        fprintf(stderr,
            "alloc failed: cap=%zu, fft_n=%zu, ~%.0f MB needed\n",
            cap, fft_n,
            (4.0 * cap * sizeof(double) + fft_n * sizeof(cplx)) / 1.0e6);
        return 1;
    }

    /* Level 0 */
    mu[0] = q;
    mu[1] = 1.0 - q;
    size_t support = 2;

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);

    int fft_active = 0;
    double max_truncated = 0.0;

    for (int k = 0; k < K; k++) {
        /* Untruncated autoconv has support [0, 2*support - 2]; we cap at cap. */
        const size_t conv_target = (2 * support - 1 < cap) ? 2 * support - 1 : cap;

        /* Right-tail cdf */
        barF[support] = 0.0;
        for (ssize_t j = (ssize_t)support - 1; j >= 0; j--)
            barF[j] = barF[j + 1] + mu[j];

        /* Pick algorithm */
        const int use_fft = !no_fft && (fft_active || support > 4096);
        if (use_fft && !fft_active) {
            fft_active = 1;
            fprintf(stderr,
                "  switching to FFT at level k=%d (support=%zu, fft_n=%zu)\n",
                k, support, fft_n);
        }

        if (use_fft) {
            real_autoconv_fft(mu, support, conv_buf, conv_target, fft_scr, fft_n);
        } else {
            real_autoconv_basic(mu, support, conv_buf, conv_target);
        }

        /* Combine: mu_new[j] = p * conv[j] + (1-p) * (barF[j]^2 - barF[j+1]^2) */
        for (size_t j = 0; j < conv_target; j++) {
            const double minp_j = (j < support)
                ? (barF[j] * barF[j] - barF[j + 1] * barF[j + 1])
                : 0.0;
            mu_new[j] = p * conv_buf[j] + (1.0 - p) * minp_j;
        }

        /* Track total mass to detect truncation. */
        double total = 0.0;
        for (size_t j = 0; j < conv_target; j++) total += mu_new[j];
        const double truncated = 1.0 - total;
        if (truncated > max_truncated) max_truncated = truncated;

        /* Swap and zero the working buffer over its used range. */
        double *tmp = mu; mu = mu_new; mu_new = tmp;
        memset(mu_new, 0, cap * sizeof(double));
        support = conv_target;

        if (K - k <= 4 || k % 5 == 0)
            fprintf(stderr,
                "  level %d -> support %zu  total=%.15g  truncated=%.3e  alg=%s\n",
                k + 1, support, total, truncated, use_fft ? "fft" : "basic");
    }

    clock_gettime(CLOCK_MONOTONIC, &t1);
    const double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    double total = 0.0;
    for (size_t j = 0; j < support; j++) total += mu[j];

    /* Write CSV */
    FILE *f = fopen(out_path, "w");
    if (!f) { perror(out_path); return 1; }
    fprintf(f, "# minplus precomputed leaf distribution\n");
    fprintf(f, "# k=%d  p=%.10f  q=%.10f  support=%zu  total=%.15g  "
               "buffer_bits=%d  max_truncated=%.3e  wall=%.3fs\n",
            K, p, q, support, total, buffer_bits, max_truncated, secs);
    fprintf(f, "value,prob\n");
    for (size_t j = 0; j < support; j++) {
        if (mu[j] > 0.0) fprintf(f, "%zu,%.18g\n", j, mu[j]);
    }
    fclose(f);

    fprintf(stderr,
        "wrote %s  k=%d  p=%.6f  q=%.6f  support=%zu  total=%.15g  "
        "max_truncated=%.3e  %.3fs\n",
        out_path, K, p, q, support, total, max_truncated, secs);

    free(mu); free(mu_new); free(barF); free(conv_buf); free(fft_scr);
    return 0;
}
