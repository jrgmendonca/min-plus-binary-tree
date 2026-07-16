# min-plus-binary-tree

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21385814.svg)](https://doi.org/10.5281/zenodo.21385814)

Monte Carlo simulator, exact leaf-table precompute, and analysis scripts for the
**min-plus process on the binary tree** — the max-type recursive distributional
equation in which each internal node of a complete binary tree applies, with
probability $p$, the operator $+$ and, with probability $1-p$, the operator
$\min$ to the values of its two children.

Companion code and data for:

> J. R. G. Mendonça, *Finite-depth scaling and an exact Bernoulli-leaf identity
> for the min-plus process on the binary tree* (2026)

## Contents

| Path | What it is |
|------|------------|
| `minplus.c` | The simulator (C11, pthreads): depth-first traversal of the upper tree over a level-wise recursion, with several leaf-distribution modes |
| `precompute_leaves.c` | Exact propagation of the full marginal µ_k through k levels, with FFT autoconvolution above support 2^12; writes cumulative leaf tables (CSV) |
| `minplus_ref.c` | Slow reference implementation used for validation |
| `Makefile` | `make` builds all three; tuned for Apple silicon — on x86-64 use `make ARCH=-march=native` |
| `precompute_all_M22.sh` | One-shot precompute of every leaf table used in the paper (27 tables, buffer 2^22; ~15 min) |
| `run_chain_M22.sh` | One-shot reproduction of the full production run set (~24 h) |
| `run_production.sh` | Named production sweeps (`precision_step4`, `noncritical_K35`, `bern_K35`, `bern_nstability`, ...) |
| `run_nearcritical_pilot.sh` | The near-critical window sweep p_N = 1/2 − c/N of Section 3.3 |
| `analyze.R`, `redo_figures.R`, `fig_nearcritical.R` | Analysis and figure generation (base R, no packages) |
| `runs/` | Histograms (`*.hist.csv`) and summaries (`*.summary.txt`) of every production run reported in the paper |
| `analysis_outputs/` | Text summaries from which the paper's tables are derived |

Per-sample streams (`*.samples.dat`) and the precomputed leaf tables (`tables/`,
~5 GB) are not shipped: both regenerate exactly — the tables because the
precompute is deterministic, the samples because every run's seed is recorded in
its `summary.txt` and in the run scripts.

## Build

```sh
make               # clang or gcc, C11; produces minplus, minplus_ref, precompute_leaves
make ARCH=-march=native   # on non-Apple hardware
```

## Quick start

Simulate the critical process (p = 1/2, all-ones leaves) at depth 20:

```sh
./minplus --depth 20 --prob 0.5 --leaves ones --samples 100000 \
          --seed 0x20000020 --out runs/critical/N20_M100000_p0.5_ones
```

Deep effective depths use a precomputed leaf table: a run at simulator depth D
with a K-level table has **effective depth N = D + K**. For example, effective
N = 60:

```sh
./precompute_leaves --prob 0.5 --q 0 --levels 35 --buffer-bits 22 \
                    --out tables/mu_K35_p0.5_q0.csv
./minplus --depth 25 --prob 0.5 --leaves table \
          --leaf-table tables/mu_K35_p0.5_q0.csv \
          --samples 100000 --seed 0x20000060 \
          --out runs/critical/N60_M100000_p0.5_K35
```

Each run writes `<prefix>.summary.txt` (metadata, moments, small-k atoms),
`<prefix>.hist.csv` (`value,count`), and `<prefix>.samples.dat` (one root value
per line; suppress with `--no-samples`).

## Reproducing the paper

```sh
./precompute_all_M22.sh          # all leaf tables (~15 min)
./run_chain_M22.sh               # all production sweeps (~24 h, N=60 cells dominate)
./run_nearcritical_pilot.sh      # near-critical window of Section 3.3 (~24 h)
Rscript redo_figures.R           # Figures 2-5, 7-8 -> figures/
Rscript fig_nearcritical.R       # Figure 6 -> figures/
```

Rough wall times on an Apple M1 Pro (10 threads, M = 1e5 samples): leaf-table
precompute at K = 35 takes ~20–60 s per table; simulation cells at D ≤ 13 run in
seconds, D = 21 in ~30 min, and D = 25 in ~7 h.

## Data formats

Leaf tables (from `precompute_leaves`): two `#` metadata lines, a `value,prob`
header, then one row per nonzero entry of mu_k, values ascending. `minplus
--leaves table` renormalizes on load, so tiny truncation deficits are harmless.

Run outputs: `*.hist.csv` is `value,count` with a `#` metadata line recording
the simulated depth, p, q, M, the leaves code, the seed, and the wall time;
`*.summary.txt` repeats the metadata together with moments and the atoms
P(X_N = k) for k ≤ 10. The `leaves=` code is the enum index: 0 ones, 1 onetwo,
2 four, 3 bern, 4 bern-shift1, 5 table.

## Reusing the pieces

- The exact one-step recursion for the full marginal (the level loop in
  `precompute_leaves.c`) is self-contained: to propagate a different integer
  leaf law, change the two lines that set mu_0.
- The simulator kernel (`draw_leaf`, `sim_levelwise`, `sim_subtree` in
  `minplus.c`, ~60 lines) depends only on the xoshiro block and `leaf_dist_t`
  and lifts as a unit; a different random operator pair replaces the single
  branch in the reduction step.
- Any distribution on the nonnegative integers can be plugged in as leaves via
  the table CSV contract above — the table need not come from
  `precompute_leaves`.
- Reproducibility: a run is determined by `--seed` and the thread count; every
  production seed is recorded in the corresponding `summary.txt` and in the
  run scripts.

## License

MIT — see `LICENSE`. If you use this code, please cite the paper (see
`CITATION.cff`).
