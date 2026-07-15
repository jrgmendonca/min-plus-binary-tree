/*
 * minplus.c - depth-first parallel simulator for the min-plus process on a
 * binary tree, and the Bernoulli-initial-values variant.
 *
 * Model: on a complete binary tree of depth D, the leaves carry i.i.d.
 * values X_0(i), and every internal node combines its two children with the
 * operator `+' with probability p and `min' with probability 1-p,
 * independently at each node. The program draws M i.i.d. samples of the
 * root value and writes their histogram, summary statistics, and
 * (optionally) the raw sample stream.
 *
 * Depth convention: --depth D is the simulated depth. With --leaves table,
 * leaves are drawn from a level-k marginal mu_k produced by
 * precompute_leaves.c, so the root is distributed as the root of a tree of
 * *effective* depth N = D + k. The built-in modes onetwo, four, and
 * bern-shift1 are the closed-form k = 1, 2, 1 instances of the same device
 * (Appendix A.1 of the companion paper).
 *
 * Memory: O(N) recursion stack plus one reusable 2^cutoff buffer per worker
 *   (--cutoff, default min(20, D), is the depth at or below which subtrees
 *   are reduced level-by-level in the buffer instead of by recursion);
 *   independent of 2^D.
 * PRNG: xoshiro256++ seeded via splitmix64 (public-domain generators by
 *   Blackman and Vigna); each worker takes its own long-jumped stream, so a
 *   run is exactly reproducible given the same --seed and thread count.
 * Output: <prefix>.hist.csv ("value,count" rows), <prefix>.summary.txt,
 *   and, unless --no-samples, <prefix>.samples.dat (one root value per line).
 *
 * Build:   make
 * Usage:   minplus --depth D --prob P --samples M [--q Q] [--threads T]
 *                  [--cutoff K] [--seed S] [--out PREFIX] [--no-samples]
 *                  --leaves {ones|onetwo|four|bern|bern-shift1|table}
 *                  [--leaf-table FILE]
 */
#include <errno.h>
#include <inttypes.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

/* ---------- xoshiro256++ ---------- */

typedef struct { uint64_t s[4]; } xoshiro_t;

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

/* xoshiro256++ "long jump" - 2^192 calls equivalent. Use to give each thread
 * its own non-overlapping stream. */
static void xoshiro_long_jump(xoshiro_t *r) {
    static const uint64_t JUMP[] = {
        0x76e15d3efefdcbbfULL, 0xc5004e441c522fb3ULL,
        0x77710069854ee241ULL, 0x39109bb02acbe635ULL };
    uint64_t s0=0,s1=0,s2=0,s3=0;
    for (int i = 0; i < 4; i++) {
        for (int b = 0; b < 64; b++) {
            if (JUMP[i] & ((uint64_t)1 << b)) {
                s0 ^= r->s[0]; s1 ^= r->s[1]; s2 ^= r->s[2]; s3 ^= r->s[3];
            }
            xoshiro_next(r);
        }
    }
    r->s[0]=s0; r->s[1]=s1; r->s[2]=s2; r->s[3]=s3;
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

/* ---------- leaf distributions ---------- */

/* The integer values of this enum (ones=0, onetwo=1, four=2, bern=3,
 * bern-shift1=4, table=5) are what the `leaves=` field means in the
 * histogram, summary, and samples file headers. */
typedef enum { LEAVES_ONES, LEAVES_ONETWO, LEAVES_FOUR, LEAVES_BERN, LEAVES_BERN_SHIFT1, LEAVES_TABLE } leaf_mode_t;

typedef struct {
    leaf_mode_t mode;
    /* Small-K modes (K <= 4): inline fixed-size arrays. */
    uint32_t values[4];
    uint64_t cum[4];   /* cumulative integer thresholds in [0, 2^64) */
    int K;
    /* LEAVES_TABLE: dynamically allocated arrays of length table_K, sorted by value
     * ascending. Sampling uses binary search through table_cum. Used to plug in a
     * precomputed level-k distribution from precompute_leaves.c, effectively shifting
     * the simulated tree depth by +k. */
    int table_K;
    uint32_t *table_values;
    uint64_t *table_cum;
} leaf_dist_t;

static uint64_t prob_to_threshold(double p) {
    if (p <= 0.0) return 0;
    if (p >= 1.0) return UINT64_MAX;
    return (uint64_t)(p * 18446744073709551616.0);
}

static void leaf_dist_init(leaf_dist_t *d, leaf_mode_t mode, double p, double q) {
    d->mode = mode;
    d->table_K = 0; d->table_values = NULL; d->table_cum = NULL;
    switch (mode) {
    case LEAVES_ONES:
        d->K = 1; d->values[0] = 1; break;
    case LEAVES_ONETWO:
        d->K = 2; d->values[0] = 1; d->values[1] = 2;
        d->cum[0] = prob_to_threshold(1.0 - p); break;
    case LEAVES_FOUR: {
        double pr1 = (1-p)*(1-p)*(1-p) + 2*p*(1-p)*(1-p);
        double pr2 = p*(1-p)*(1-p) + (1-p)*p*p;
        double pr3 = 2*(1-p)*p*p;
        d->K = 4;
        d->values[0] = 1; d->values[1] = 2; d->values[2] = 3; d->values[3] = 4;
        d->cum[0] = prob_to_threshold(pr1);
        d->cum[1] = prob_to_threshold(pr1 + pr2);
        d->cum[2] = prob_to_threshold(pr1 + pr2 + pr3);
        break;
    }
    case LEAVES_BERN:
        d->K = 2; d->values[0] = 0; d->values[1] = 1;
        d->cum[0] = prob_to_threshold(q); break;
    case LEAVES_BERN_SHIFT1: {
        double pr0 = q*q + 2*q*(1-q)*(1-p);
        double pr1 = 2*q*(1-q)*p + (1-q)*(1-q)*(1-p);
        d->K = 3;
        d->values[0] = 0; d->values[1] = 1; d->values[2] = 2;
        d->cum[0] = prob_to_threshold(pr0);
        d->cum[1] = prob_to_threshold(pr0 + pr1);
        break;
    }
    case LEAVES_TABLE:
        /* Caller must use leaf_dist_init_from_file instead. */
        d->K = 0;
        break;
    }
}

/* Load a precomputed mu_k distribution from the CSV produced by
 * precompute_leaves.c. Skips comment lines (#...) and the header row
 * "value,prob". Each remaining row is "<value>,<prob>". The values do not
 * have to be contiguous and may have prob == 0 entries omitted. */
static int leaf_dist_init_from_file(leaf_dist_t *d, const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return -1; }
    /* Pass 1: count data rows. */
    char line[1024];
    int n = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (!strncmp(line, "value,", 6)) continue;
        n++;
    }
    if (n <= 0) {
        fprintf(stderr, "leaf table %s: no data rows\n", path);
        fclose(f); return -1;
    }
    /* Pass 2: parse. */
    rewind(f);
    uint32_t *vals  = (uint32_t *)malloc(n * sizeof(uint32_t));
    double   *probs = (double   *)malloc(n * sizeof(double));
    if (!vals || !probs) { fprintf(stderr, "alloc\n"); fclose(f); return -1; }
    int i = 0;
    while (fgets(line, sizeof line, f) && i < n) {
        if (line[0] == '#' || line[0] == '\n') continue;
        if (!strncmp(line, "value,", 6)) continue;
        unsigned long long v;
        double pr;
        if (sscanf(line, "%llu,%lf", &v, &pr) != 2) {
            fprintf(stderr, "leaf table %s: malformed line: %s", path, line);
            free(vals); free(probs); fclose(f); return -1;
        }
        vals[i]  = (uint32_t)v;
        probs[i] = pr;
        i++;
    }
    fclose(f);

    /* Verify ordering and normalize. The CSV from precompute_leaves is sorted
     * by value, but be defensive: insertion-sort for safety. */
    for (int k = 1; k < n; k++) {
        if (vals[k] >= vals[k-1]) continue;
        uint32_t v = vals[k]; double pr = probs[k];
        int j = k - 1;
        while (j >= 0 && vals[j] > v) { vals[j+1] = vals[j]; probs[j+1] = probs[j]; j--; }
        vals[j+1] = v; probs[j+1] = pr;
    }
    double total = 0.0;
    for (int k = 0; k < n; k++) total += probs[k];
    if (total <= 0.0) {
        fprintf(stderr, "leaf table %s: total mass = %g\n", path, total);
        free(vals); free(probs); return -1;
    }

    /* Build cumulative integer thresholds. cum[k] is the upper-exclusive
     * bound for outcome vals[k] in [0, 2^64). cum[n-1] = UINT64_MAX. */
    d->mode      = LEAVES_TABLE;
    d->K         = 0;             /* 0 means "use table_*" in draw_leaf */
    d->table_K   = n;
    d->table_values = vals;
    d->table_cum    = (uint64_t *)malloc(n * sizeof(uint64_t));
    if (!d->table_cum) { fprintf(stderr, "alloc\n"); free(vals); free(probs); return -1; }
    double acc = 0.0;
    for (int k = 0; k < n - 1; k++) {
        acc += probs[k];
        d->table_cum[k] = prob_to_threshold(acc / total);
    }
    d->table_cum[n - 1] = UINT64_MAX;
    free(probs);
    fprintf(stderr, "loaded leaf table %s  (entries=%d  total=%.15g)\n",
            path, n, total);
    return 0;
}

static void leaf_dist_free(leaf_dist_t *d) {
    free(d->table_values); d->table_values = NULL;
    free(d->table_cum);    d->table_cum    = NULL;
    d->table_K = 0;
}

static inline uint32_t draw_leaf(xoshiro_t *r, const leaf_dist_t *d) {
    if (d->K == 1) return d->values[0];
    uint64_t u = xoshiro_next(r);
    if (d->mode == LEAVES_TABLE) {
        /* Binary search through the cumulative integer thresholds. */
        int lo = 0, hi = d->table_K - 1;
        while (lo < hi) {
            int mid = (lo + hi) >> 1;
            if (u < d->table_cum[mid]) hi = mid;
            else lo = mid + 1;
        }
        return d->table_values[lo];
    }
    if (d->K == 2) return d->values[u < d->cum[0] ? 0 : 1];
    /* K=3 or K=4 */
    if (u < d->cum[0]) return d->values[0];
    if (u < d->cum[1]) return d->values[1];
    if (d->K == 3 || u < d->cum[2]) return d->values[2];
    return d->values[3];
}

/* ---------- hybrid simulator ----------
 *
 * The bottom K levels of a subtree are reduced level-by-level on a reusable
 * per-thread buffer of size 2^K (tight inner loop, good for the M1's wide
 * cores). Above K we recurse, so total memory is 2^K per thread independent
 * of N. Default K = min(20, N) -> 8 MB per worker.
 */

__attribute__((always_inline))
static inline uint64_t sim_levelwise(int K, uint64_t p_thresh, const leaf_dist_t *d,
                                     xoshiro_t * __restrict__ r,
                                     uint64_t * __restrict__ buf) {
    size_t k = (size_t)1 << K;
    if (d->K == 1) {
        uint64_t v = d->values[0];
        for (size_t i = 0; i < k; i++) buf[i] = v;
    } else {
        for (size_t i = 0; i < k; i++) buf[i] = draw_leaf(r, d);
    }
    while (k > 1) {
        k >>= 1;
        for (size_t i = 0; i < k; i++) {
            uint64_t a = buf[2*i], b = buf[2*i + 1];
            buf[i] = (xoshiro_next(r) < p_thresh) ? (a + b) : (a < b ? a : b);
        }
    }
    return buf[0];
}

static uint64_t sim_subtree(int depth, int K, uint64_t p_thresh, const leaf_dist_t *d,
                            xoshiro_t * __restrict__ r, uint64_t * __restrict__ buf) {
    if (depth <= K) return sim_levelwise(depth, p_thresh, d, r, buf);
    uint64_t a = sim_subtree(depth - 1, K, p_thresh, d, r, buf);
    uint64_t b = sim_subtree(depth - 1, K, p_thresh, d, r, buf);
    return (xoshiro_next(r) < p_thresh) ? (a + b) : (a < b ? a : b);
}

/* ---------- histogram ---------- */

#define HIST_FAST_MAX (1 << 16)   /* values 0..65535 in flat array */

typedef struct {
    uint64_t fast[HIST_FAST_MAX];
    /* Sparse overflow: parallel arrays of (value, count), sorted by value. */
    uint64_t *over_val;
    uint64_t *over_cnt;
    size_t over_n, over_cap;
} hist_t;

static void hist_init(hist_t *h) {
    memset(h->fast, 0, sizeof h->fast);
    h->over_val = NULL; h->over_cnt = NULL; h->over_n = 0; h->over_cap = 0;
}

static void hist_free(hist_t *h) {
    free(h->over_val); free(h->over_cnt);
}

static void hist_add_overflow(hist_t *h, uint64_t v) {
    /* Linear probe; overflow values should be very rare. */
    for (size_t i = 0; i < h->over_n; i++) {
        if (h->over_val[i] == v) { h->over_cnt[i]++; return; }
    }
    if (h->over_n == h->over_cap) {
        size_t nc = h->over_cap ? h->over_cap * 2 : 64;
        h->over_val = (uint64_t *)realloc(h->over_val, nc * sizeof(uint64_t));
        h->over_cnt = (uint64_t *)realloc(h->over_cnt, nc * sizeof(uint64_t));
        if (!h->over_val || !h->over_cnt) { fprintf(stderr, "hist realloc failed\n"); exit(1); }
        h->over_cap = nc;
    }
    h->over_val[h->over_n] = v; h->over_cnt[h->over_n] = 1; h->over_n++;
}

static inline void hist_add(hist_t *h, uint64_t v) {
    if (v < HIST_FAST_MAX) h->fast[v]++;
    else hist_add_overflow(h, v);
}

static void hist_merge(hist_t *dst, const hist_t *src) {
    for (size_t i = 0; i < HIST_FAST_MAX; i++) dst->fast[i] += src->fast[i];
    for (size_t i = 0; i < src->over_n; i++) {
        uint64_t v = src->over_val[i], c = src->over_cnt[i];
        size_t j;
        for (j = 0; j < dst->over_n; j++) if (dst->over_val[j] == v) { dst->over_cnt[j] += c; break; }
        if (j < dst->over_n) continue;
        if (dst->over_n == dst->over_cap) {
            size_t nc = dst->over_cap ? dst->over_cap * 2 : 64;
            dst->over_val = (uint64_t *)realloc(dst->over_val, nc * sizeof(uint64_t));
            dst->over_cnt = (uint64_t *)realloc(dst->over_cnt, nc * sizeof(uint64_t));
            if (!dst->over_val || !dst->over_cnt) { fprintf(stderr, "merge realloc failed\n"); exit(1); }
            dst->over_cap = nc;
        }
        dst->over_val[dst->over_n] = v; dst->over_cnt[dst->over_n] = c; dst->over_n++;
    }
}

static int cmp_u64(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

/* ---------- worker ---------- */

typedef struct {
    int depth;
    int cutoff_K;
    uint64_t p_thresh;
    const leaf_dist_t *leaves;
    uint64_t n_samples;
    xoshiro_t rng;
    uint64_t *buf;   /* reusable 2^K buffer for level-wise reduction */
    /* outputs */
    hist_t hist;
    uint64_t *samples;     /* may be NULL if not collecting raw samples */
    uint64_t sample_off;   /* offset into shared samples array */
    /* stats */
    long double sum, sumsq;
    uint64_t mn, mx;
} worker_t;

static void *worker_run(void *arg) {
    worker_t *w = (worker_t *)arg;
    hist_init(&w->hist);
    w->sum = 0; w->sumsq = 0; w->mn = UINT64_MAX; w->mx = 0;
    for (uint64_t s = 0; s < w->n_samples; s++) {
        uint64_t v = sim_subtree(w->depth, w->cutoff_K, w->p_thresh, w->leaves, &w->rng, w->buf);
        hist_add(&w->hist, v);
        w->sum += v; w->sumsq += (long double)v * v;
        if (v < w->mn) w->mn = v; if (v > w->mx) w->mx = v;
        if (w->samples) w->samples[w->sample_off + s] = v;
    }
    return NULL;
}

/* sqrt without pulling in math.h; Newton's method, fine for telemetry. */
static double sqrt_(double x) {
    if (x <= 0) return 0;
    double r = x;
    for (int i = 0; i < 32; i++) r = 0.5 * (r + x / r);
    return r;
}

/* ---------- output ---------- */

static void write_histogram(const char *path, const hist_t *h, int N, double p, double q,
                            uint64_t M, leaf_mode_t mode, uint64_t seed, double secs) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    fprintf(f, "# minplus histogram\n");
    fprintf(f, "# N=%d  p=%.10f  q=%.10f  M=%" PRIu64 "  leaves=%d  seed=0x%016" PRIx64 "  wall=%.3fs\n",
            N, p, q, M, (int)mode, seed, secs);
    fprintf(f, "value,count\n");
    for (size_t v = 0; v < HIST_FAST_MAX; v++)
        if (h->fast[v]) fprintf(f, "%zu,%" PRIu64 "\n", v, h->fast[v]);
    /* Sort and emit overflow */
    if (h->over_n) {
        /* Build sorted index of overflow */
        size_t *idx = (size_t *)malloc(h->over_n * sizeof(size_t));
        for (size_t i = 0; i < h->over_n; i++) idx[i] = i;
        /* Simple sort by value */
        for (size_t i = 1; i < h->over_n; i++) {
            size_t j = i;
            while (j > 0 && h->over_val[idx[j-1]] > h->over_val[idx[j]]) {
                size_t t = idx[j-1]; idx[j-1] = idx[j]; idx[j] = t; j--;
            }
        }
        for (size_t i = 0; i < h->over_n; i++)
            fprintf(f, "%" PRIu64 ",%" PRIu64 "\n", h->over_val[idx[i]], h->over_cnt[idx[i]]);
        free(idx);
    }
    fclose(f);
}

static void write_summary(const char *path, int N, double p, double q, uint64_t M,
                          leaf_mode_t mode, uint64_t seed, int threads, double secs,
                          double mean, double var, uint64_t mn, uint64_t mx,
                          const hist_t *h) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    fprintf(f, "minplus simulation summary\n");
    fprintf(f, "  depth N        : %d\n", N);
    fprintf(f, "  prob p         : %.10f\n", p);
    fprintf(f, "  prob q         : %.10f\n", q);
    fprintf(f, "  samples M      : %" PRIu64 "\n", M);
    fprintf(f, "  leaves         : %d\n", (int)mode);
    fprintf(f, "  threads        : %d\n", threads);
    fprintf(f, "  seed           : 0x%016" PRIx64 "\n", seed);
    fprintf(f, "  wall time      : %.3fs  (%.0f samples/s)\n", secs, (double)M / secs);
    fprintf(f, "  mean(X_N)      : %.6f\n", mean);
    fprintf(f, "  var(X_N)       : %.6f\n", var);
    fprintf(f, "  std(X_N)       : %.6f\n", var > 0 ? sqrt_(var) : 0.0);
    fprintf(f, "  min            : %" PRIu64 "\n", mn);
    fprintf(f, "  max            : %" PRIu64 "\n", mx);
    fprintf(f, "  P(X_N=k) for k=0..10:\n");
    for (int k = 0; k <= 10; k++) {
        uint64_t c = (k < HIST_FAST_MAX) ? h->fast[k] : 0;
        fprintf(f, "    k=%-3d  count=%-12" PRIu64 "  prob=%.6e\n", k, c, (double)c / M);
    }
    fclose(f);
}

/* ---------- main ---------- */

static void usage(const char *p) {
    fprintf(stderr,
        "usage: %s --depth N --prob P --samples M [--q Q] [--threads T] [--cutoff K]\n"
        "          [--seed S] --leaves {ones|onetwo|four|bern|bern-shift1|table}\n"
        "          [--leaf-table FILE] [--out PREFIX] [--no-samples]\n"
        "  --leaves table requires --leaf-table FILE (CSV from precompute_leaves).\n", p);
    exit(2);
}

int main(int argc, char **argv) {
    int N = -1, threads = 0, want_samples = 1, cutoff_K = -1;
    uint64_t M = 0;
    double p = -1.0, q = 0.0;
    uint64_t seed = 0; int seed_set = 0;
    leaf_mode_t mode = LEAVES_ONES; int mode_set = 0;
    const char *out_prefix = "minplus_run";
    const char *leaf_table_path = NULL;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (!strcmp(a, "--depth") && i+1 < argc) N = atoi(argv[++i]);
        else if (!strcmp(a, "--prob") && i+1 < argc) p = atof(argv[++i]);
        else if (!strcmp(a, "--q") && i+1 < argc) q = atof(argv[++i]);
        else if (!strcmp(a, "--samples") && i+1 < argc) M = strtoull(argv[++i], NULL, 0);
        else if (!strcmp(a, "--threads") && i+1 < argc) threads = atoi(argv[++i]);
        else if (!strcmp(a, "--cutoff") && i+1 < argc) cutoff_K = atoi(argv[++i]);
        else if (!strcmp(a, "--seed") && i+1 < argc) { seed = strtoull(argv[++i], NULL, 0); seed_set = 1; }
        else if (!strcmp(a, "--out") && i+1 < argc) out_prefix = argv[++i];
        else if (!strcmp(a, "--no-samples")) want_samples = 0;
        else if (!strcmp(a, "--leaf-table") && i+1 < argc) leaf_table_path = argv[++i];
        else if (!strcmp(a, "--leaves") && i+1 < argc) {
            const char *m = argv[++i]; mode_set = 1;
            if (!strcmp(m, "ones")) mode = LEAVES_ONES;
            else if (!strcmp(m, "onetwo")) mode = LEAVES_ONETWO;
            else if (!strcmp(m, "four")) mode = LEAVES_FOUR;
            else if (!strcmp(m, "bern")) mode = LEAVES_BERN;
            else if (!strcmp(m, "bern-shift1")) mode = LEAVES_BERN_SHIFT1;
            else if (!strcmp(m, "table")) mode = LEAVES_TABLE;
            else { fprintf(stderr, "unknown leaves mode: %s\n", m); usage(argv[0]); }
        }
        else usage(argv[0]);
    }
    if (N < 1 || M < 1 || p < 0.0 || p > 1.0 || !mode_set) usage(argv[0]);
    if (mode == LEAVES_TABLE && !leaf_table_path) {
        fprintf(stderr, "--leaves table requires --leaf-table FILE\n"); usage(argv[0]);
    }
    if (!seed_set) seed = (uint64_t)time(NULL) ^ ((uint64_t)getpid() << 32);
    if (threads <= 0) {
        long n = sysconf(_SC_NPROCESSORS_ONLN);
        threads = (n > 0) ? (int)n : 4;
        if (threads > 16) threads = 16;
    }
    if ((uint64_t)threads > M) threads = (int)M;
    if (cutoff_K < 0) cutoff_K = (N < 20) ? N : 20;
    if (cutoff_K > N) cutoff_K = N;

    leaf_dist_t leaves;
    leaf_dist_init(&leaves, mode, p, q);
    if (mode == LEAVES_TABLE) {
        if (leaf_dist_init_from_file(&leaves, leaf_table_path) != 0) return 1;
    }
    uint64_t p_thresh = prob_to_threshold(p);

    /* Per-thread samples partition. */
    worker_t *ws = (worker_t *)calloc(threads, sizeof(worker_t));
    pthread_t *ths = (pthread_t *)calloc(threads, sizeof(pthread_t));
    uint64_t *all_samples = NULL;
    if (want_samples) {
        all_samples = (uint64_t *)malloc(M * sizeof(uint64_t));
        if (!all_samples) { fprintf(stderr, "samples malloc failed (%" PRIu64 ")\n", M); return 1; }
    }

    uint64_t per = M / threads, rem = M % threads, off = 0;
    xoshiro_t base; xoshiro_seed(&base, seed);
    size_t buf_sz = (size_t)1 << cutoff_K;
    for (int t = 0; t < threads; t++) {
        ws[t].depth = N;
        ws[t].cutoff_K = cutoff_K;
        ws[t].p_thresh = p_thresh;
        ws[t].leaves = &leaves;
        ws[t].n_samples = per + (t < (int)rem ? 1 : 0);
        ws[t].rng = base;
        for (int j = 0; j < t; j++) xoshiro_long_jump(&ws[t].rng);  /* unique stream per worker */
        ws[t].buf = (uint64_t *)malloc(buf_sz * sizeof(uint64_t));
        if (!ws[t].buf) { fprintf(stderr, "thread buf malloc failed (K=%d)\n", cutoff_K); return 1; }
        ws[t].samples = all_samples;
        ws[t].sample_off = off;
        off += ws[t].n_samples;
    }

    struct timespec t0, t1;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int t = 0; t < threads; t++) {
        if (pthread_create(&ths[t], NULL, worker_run, &ws[t]) != 0) {
            fprintf(stderr, "pthread_create failed: %s\n", strerror(errno)); return 1;
        }
    }
    for (int t = 0; t < threads; t++) pthread_join(ths[t], NULL);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    double secs = (t1.tv_sec - t0.tv_sec) + (t1.tv_nsec - t0.tv_nsec) / 1e9;

    /* Merge histograms and stats. */
    hist_t H; hist_init(&H);
    long double sum = 0, sumsq = 0; uint64_t mn = UINT64_MAX, mx = 0;
    for (int t = 0; t < threads; t++) {
        hist_merge(&H, &ws[t].hist);
        sum += ws[t].sum; sumsq += ws[t].sumsq;
        if (ws[t].mn < mn) mn = ws[t].mn;
        if (ws[t].mx > mx) mx = ws[t].mx;
    }
    double mean = (double)(sum / M);
    double var = (double)(sumsq / M) - mean * mean;

    fprintf(stderr,
        "# minplus  N=%d  p=%.6f  q=%.6f  M=%" PRIu64 "  leaves=%d  threads=%d  seed=0x%016" PRIx64 "\n"
        "# wall=%.3fs  rate=%.0f samples/s  mean=%.6f  var=%.6f  min=%" PRIu64 "  max=%" PRIu64 "\n",
        N, p, q, M, (int)mode, threads, seed, secs, (double)M / secs, mean, var, mn, mx);

    char path[1024];
    snprintf(path, sizeof path, "%s.hist.csv", out_prefix);
    write_histogram(path, &H, N, p, q, M, mode, seed, secs);
    snprintf(path, sizeof path, "%s.summary.txt", out_prefix);
    write_summary(path, N, p, q, M, mode, seed, threads, secs, mean, var, mn, mx, &H);
    if (want_samples) {
        snprintf(path, sizeof path, "%s.samples.dat", out_prefix);
        FILE *f = fopen(path, "w");
        if (!f) { perror(path); return 1; }
        fprintf(f, "# N=%d p=%.6f q=%.6f M=%" PRIu64 " leaves=%d seed=0x%016" PRIx64 "\n",
                N, p, q, M, (int)mode, seed);
        for (uint64_t s = 0; s < M; s++) fprintf(f, "%" PRIu64 "\n", all_samples[s]);
        fclose(f);
    }

    /* cleanup */
    for (int t = 0; t < threads; t++) { hist_free(&ws[t].hist); free(ws[t].buf); }
    hist_free(&H);
    leaf_dist_free(&leaves);
    free(ws); free(ths); free(all_samples);
    (void)cmp_u64;
    return 0;
}
