/*
 * minplus_ref.c - single-threaded reference implementation for
 * cross-validation (a straight port of the original R prototype).
 *
 * Direct level-by-level reduction: allocate 2^N leaves, draw one uniform per
 * internal node, apply (sum | min) by Bernoulli(p). Memory O(2^N). Use only
 * at small N (<= 26), to confirm that the optimised depth-first simulator
 * minplus.c agrees in distribution. For production runs use minplus.c.
 *
 * Build:   make minplus_ref
 * Usage:   minplus_ref --depth N --prob P [--q Q] --samples M [--seed S]
 *                      --leaves {ones|onetwo|four|bern|bern-shift1}
 *                      [--out PREFIX]
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {
    uint64_t s[4];
} xoshiro_t;

static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

static inline uint64_t xoshiro_next(xoshiro_t *r) {
    uint64_t res = rotl(r->s[0] + r->s[3], 23) + r->s[0];
    uint64_t t = r->s[1] << 17;
    r->s[2] ^= r->s[0]; r->s[3] ^= r->s[1];
    r->s[1] ^= r->s[2]; r->s[0] ^= r->s[3];
    r->s[2] ^= t;
    r->s[3] = rotl(r->s[3], 45);
    return res;
}

static inline uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static void xoshiro_seed(xoshiro_t *r, uint64_t seed) {
    uint64_t s = seed;
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&s);
}

/* Convert bottom 53 bits to a [0,1) double, then compare. We use an integer
 * threshold instead, to keep the hot path branch-only on integer compares. */
static inline int bernoulli(xoshiro_t *r, uint64_t threshold) {
    return xoshiro_next(r) < threshold;
}

/* Sample an index in [0,K) given cumulative thresholds (each in [0, 2^64)). */
static inline uint32_t sample_discrete(xoshiro_t *r, const uint64_t *cum, int K) {
    uint64_t u = xoshiro_next(r);
    for (int i = 0; i < K - 1; i++) if (u < cum[i]) return i;
    return K - 1;
}

typedef enum { LEAVES_ONES, LEAVES_ONETWO, LEAVES_FOUR, LEAVES_BERN, LEAVES_BERN_SHIFT1 } leaf_mode_t;

typedef struct {
    leaf_mode_t mode;
    uint32_t values[4];
    uint64_t cum[4];   /* cumulative integer thresholds */
    int K;
} leaf_dist_t;

static double prob_to_double(uint64_t t) { return (double)t / 18446744073709551616.0; }

static uint64_t prob_to_threshold(double p) {
    if (p <= 0.0) return 0;
    if (p >= 1.0) return UINT64_MAX;
    /* 2^64 * p, clamped */
    return (uint64_t)(p * 18446744073709551616.0);
}

static void leaf_dist_init(leaf_dist_t *d, leaf_mode_t mode, double p, double q) {
    d->mode = mode;
    switch (mode) {
    case LEAVES_ONES:
        d->K = 1; d->values[0] = 1;
        break;
    case LEAVES_ONETWO: {
        /* One-level leaf construction: X_0 in {1,2} with probs (1-p, p) --
         * effective depth +1. */
        d->K = 2; d->values[0] = 1; d->values[1] = 2;
        d->cum[0] = prob_to_threshold(1.0 - p);
        break;
    }
    case LEAVES_FOUR: {
        /* Two-level leaf construction -- effective depth +2. */
        double pr1 = (1-p)*(1-p)*(1-p) + 2*p*(1-p)*(1-p);
        double pr2 = p*(1-p)*(1-p) + (1-p)*p*p;
        double pr3 = 2*(1-p)*p*p;
        /* pr4 = p^3, implicit */
        d->K = 4;
        d->values[0] = 1; d->values[1] = 2; d->values[2] = 3; d->values[3] = 4;
        d->cum[0] = prob_to_threshold(pr1);
        d->cum[1] = prob_to_threshold(pr1 + pr2);
        d->cum[2] = prob_to_threshold(pr1 + pr2 + pr3);
        break;
    }
    case LEAVES_BERN: {
        /* Bernoulli variant: P(X_0=0)=q, P(X_0=1)=1-q */
        d->K = 2; d->values[0] = 0; d->values[1] = 1;
        d->cum[0] = prob_to_threshold(q);
        break;
    }
    case LEAVES_BERN_SHIFT1: {
        /* One-step law of Bernoulli leaves -- effective depth +1 of the
         * variant. */
        double pr0 = q*q + 2*q*(1-q)*(1-p);
        double pr1 = 2*q*(1-q)*p + (1-q)*(1-q)*(1-p);
        /* pr2 = (1-q)^2 * p, implicit */
        d->K = 3;
        d->values[0] = 0; d->values[1] = 1; d->values[2] = 2;
        d->cum[0] = prob_to_threshold(pr0);
        d->cum[1] = prob_to_threshold(pr0 + pr1);
        break;
    }
    }
}

static inline uint32_t draw_leaf(xoshiro_t *r, const leaf_dist_t *d) {
    if (d->K == 1) return d->values[0];
    return d->values[sample_discrete(r, d->cum, d->K)];
}

/* One MC sample on a tree of depth N, level-by-level (mirrors R algorithm). */
static uint64_t simulate_one(int N, uint64_t p_thresh, const leaf_dist_t *d, xoshiro_t *r) {
    size_t k = (size_t)1 << N;
    uint64_t *x = (uint64_t *)malloc(k * sizeof(uint64_t));
    if (!x) { fprintf(stderr, "malloc(%zu) failed\n", k * sizeof(uint64_t)); exit(1); }

    for (size_t i = 0; i < k; i++) x[i] = draw_leaf(r, d);

    while (k > 1) {
        k >>= 1;
        for (size_t i = 0; i < k; i++) {
            uint64_t a = x[2*i], b = x[2*i + 1];
            x[i] = bernoulli(r, p_thresh) ? (a + b) : (a < b ? a : b);
        }
    }
    uint64_t result = x[0];
    free(x);
    return result;
}

static void usage(const char *p) {
    fprintf(stderr,
        "usage: %s --depth N --prob P --samples M [--q Q] [--seed S]\n"
        "          --leaves {ones|onetwo|four|bern|bern-shift1} [--out PREFIX]\n", p);
    exit(2);
}

int main(int argc, char **argv) {
    int N = -1, M = -1;
    double p = -1.0, q = 0.0;
    uint64_t seed = 0;
    int seed_set = 0;
    leaf_mode_t mode = LEAVES_ONES;
    int mode_set = 0;
    const char *out_prefix = NULL;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--depth") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(a, "--prob") && i+1 < argc) p = atof(argv[++i]);
        else if (!strcmp(a, "--q") && i+1 < argc) q = atof(argv[++i]);
        else if (!strcmp(a, "--samples") && i+1 < argc) M = atoi(argv[++i]);
        else if (!strcmp(a, "--seed") && i+1 < argc) { seed = strtoull(argv[++i], NULL, 0); seed_set = 1; }
        else if (!strcmp(a, "--out") && i+1 < argc) out_prefix = argv[++i];
        else if (!strcmp(a, "--leaves") && i+1 < argc) {
            const char *m = argv[++i]; mode_set = 1;
            if (!strcmp(m, "ones")) mode = LEAVES_ONES;
            else if (!strcmp(m, "onetwo")) mode = LEAVES_ONETWO;
            else if (!strcmp(m, "four")) mode = LEAVES_FOUR;
            else if (!strcmp(m, "bern")) mode = LEAVES_BERN;
            else if (!strcmp(m, "bern-shift1")) mode = LEAVES_BERN_SHIFT1;
            else { fprintf(stderr, "unknown leaves mode: %s\n", m); usage(argv[0]); }
        }
        else usage(argv[0]);
    }
    if (N < 1 || M < 1 || p < 0.0 || p > 1.0 || !mode_set) usage(argv[0]);
    if (!seed_set) seed = (uint64_t)time(NULL) ^ ((uint64_t)getpid() << 32);

    leaf_dist_t leaves;
    leaf_dist_init(&leaves, mode, p, q);
    uint64_t p_thresh = prob_to_threshold(p);

    xoshiro_t rng;
    xoshiro_seed(&rng, seed);

    uint64_t *samples = (uint64_t *)malloc((size_t)M * sizeof(uint64_t));
    if (!samples) { fprintf(stderr, "malloc samples failed\n"); return 1; }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int s = 0; s < M; s++) samples[s] = simulate_one(N, p_thresh, &leaves, &rng);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    /* Summary */
    long double sum = 0, sumsq = 0; uint64_t mn = UINT64_MAX, mx = 0;
    uint64_t eq1 = 0;
    for (int s = 0; s < M; s++) {
        uint64_t v = samples[s];
        sum += v; sumsq += (long double)v * v;
        if (v < mn) mn = v; if (v > mx) mx = v;
        if (v == 1) eq1++;
    }
    double mean = (double)(sum / M);
    double var = (double)(sumsq / M) - mean*mean;

    fprintf(stderr,
        "# minplus_ref  N=%d  p=%.6f  q=%.6f  M=%d  leaves=%d  seed=0x%016" PRIx64 "\n"
        "# wall=%.3fs  mean=%.6f  var=%.6f  min=%" PRIu64 "  max=%" PRIu64 "  P(X=1)=%.6f\n",
        N, p, q, M, (int)mode, seed, secs, mean, var, mn, mx, (double)eq1 / M);

    if (out_prefix) {
        char path[1024];
        snprintf(path, sizeof path, "%s.samples.dat", out_prefix);
        FILE *f = fopen(path, "w");
        if (!f) { perror(path); return 1; }
        fprintf(f, "# N=%d p=%.6f q=%.6f M=%d leaves=%d seed=0x%016" PRIx64 "\n", N, p, q, M, (int)mode, seed);
        for (int s = 0; s < M; s++) fprintf(f, "%" PRIu64 "\n", samples[s]);
        fclose(f);
    } else {
        for (int s = 0; s < M; s++) printf("%" PRIu64 "\n", samples[s]);
    }

    free(samples);
    (void)prob_to_double;
    return 0;
}
