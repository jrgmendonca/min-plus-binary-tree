#!/usr/bin/env bash
#
# Production-run orchestrator for the min-plus simulator.
#
# Usage:
#   ./run_production.sh smoke         # ~12 min — critical N=22,24,26 (eff 24,26,28)
#   ./run_production.sh medium        # ~2.7 h  — critical N=28,30 (eff 30,32)
#   ./run_production.sh long          # ~9 h    — critical N=32 (eff 34)
#   ./run_production.sh xlong         # DISABLED — superseded by precision_step4
#   ./run_production.sh precision_step4 # ~2.5 h — eff 20, 36, 40, 44, 48, 52, 56, 60 at M=1e5
#   ./run_production.sh bern_smoke    # ~10 min — Bernoulli sweep at N=24
#   ./run_production.sh bern_n26      # ~35 min — Bernoulli sweep at N=26, q in [0,0.5]
#   ./run_production.sh bern_n26_upper # ~35 min — Bernoulli sweep at N=26, q in [0.55,1.0]
#   ./run_production.sh bern          # ~85 min — Bernoulli sweep at N=28
#   ./run_production.sh noncritical   # ~2.2 h  — p in {0.1,0.2,0.3,0.4} at N=28 (eff 30)
#   ./run_production.sh all_short     # smoke + bern_smoke + noncritical
#
# All runs use --leaves four (Eq. 22) for critical/non-critical, which shifts
# the effective tree depth by +2 relative to --depth. Bernoulli runs use the
# bare leaves at depth N (no shift). Outputs go to runs/<class>/<name>.{samples.dat,
# hist.csv,summary.txt}; logs to logs/<name>.log. Reruns skip already-completed
# runs (idempotent on summary.txt).

set -euo pipefail
cd "$(dirname "$0")"

BIN=./minplus
test -x "$BIN" || { echo "build first: make"; exit 1; }

run() {
    # run <out_prefix> <minplus args...>
    local out="$1"; shift
    local logf="logs/$(basename "$out").log"
    if [[ -f "${out}.summary.txt" ]]; then
        echo "[skip] ${out}"
        return 0
    fi
    echo "[run ] ${out}"
    local t0 t1
    t0=$(date +%s)
    "$BIN" "$@" --out "$out" >"$logf" 2>&1
    t1=$(date +%s)
    printf "[done] %s  (%ds)\n" "$out" "$((t1 - t0))"
    tail -2 "$logf" | sed 's/^/       /'
}

# --- critical p=0.5 ladder (leaves=four, effective depth = N+2) -------------

smoke() {
    run runs/critical/N22_M10000_p0.5_four \
        --depth 22 --prob 0.5 --leaves four --samples 10000 --seed 0x10000022
    run runs/critical/N24_M10000_p0.5_four \
        --depth 24 --prob 0.5 --leaves four --samples 10000 --seed 0x10000024
    run runs/critical/N26_M10000_p0.5_four \
        --depth 26 --prob 0.5 --leaves four --samples 10000 --seed 0x10000026
}

medium() {
    run runs/critical/N28_M10000_p0.5_four \
        --depth 28 --prob 0.5 --leaves four --samples 10000 --seed 0x10000028
    run runs/critical/N30_M10000_p0.5_four \
        --depth 30 --prob 0.5 --leaves four --samples 10000 --seed 0x10000030
}

long_run() {
    run runs/critical/N32_M10000_p0.5_four \
        --depth 32 --prob 0.5 --leaves four --samples 10000 --seed 0x10000032
}

xlong() {
    # N=34 four ~ effective 36. M=2500 keeps wall ~9h; combine seeds later if
    # more samples needed.
    # Disabled in favour of precision_step4, which covers eff 36 at M=1e5
    # via the mu_K=35 leaf table for ~50x better statistics in less time.
    run runs/critical/N34_M2500_p0.5_four \
        --depth 34 --prob 0.5 --leaves four --samples 2500  --seed 0x10000034
}

precision_step4() {
    # High-precision step-by-4 grid at M=1e5 using two precomputed leaf tables:
    # mu_K=20 for eff 24, 28, 32 (D = N - 20 in {4, 8, 12}, all sub-second);
    # mu_K=35 for eff 36, 40, 44, 48, 52, 56, 60 (D = N - 35 in {1, ..., 25},
    # dominated by eff 60 at ~140 min). Plus eff 20 baseline at depth 20 with
    # no leaf shift. Wall budget: ~2.5 h end-to-end.
    local M=100000
    local TABLE_LO=tables/mu_K20_p0.5_q0.csv
    local TABLE_HI=tables/mu_K35_p0.5_q0.csv
    for t in "$TABLE_LO" "$TABLE_HI"; do
        if [[ ! -f "$t" ]]; then
            echo "missing leaf table $t; run precompute_leaves first" >&2
            exit 1
        fi
    done
    # eff 20: no shift, all-ones leaves.
    run "runs/critical/N20_M${M}_p0.5_ones" \
        --depth 20 --prob 0.5 --leaves ones --samples "$M" --seed "0x20000020"
    # eff 24, 28, 32 via mu_20 (D = N - 20).
    local D effN
    for D in 4 8 12; do
        effN=$((D + 20))
        run "runs/critical/N${effN}_M${M}_p0.5_K20" \
            --depth "$D" --prob 0.5 --leaves table --leaf-table "$TABLE_LO" \
            --samples "$M" --seed "0x20000${effN}"
    done
    # eff 36, 40, 44, ..., 60 via mu_35 (D = N - 35).
    for D in 1 5 9 13 17 21 25; do
        effN=$((D + 35))
        run "runs/critical/N${effN}_M${M}_p0.5_K35" \
            --depth "$D" --prob 0.5 --leaves table --leaf-table "$TABLE_HI" \
            --samples "$M" --seed "0x20000${effN}"
    done
}

# --- Bernoulli variant (no shift; bare bern leaves) -------------------------

bern_smoke() {
    # Quick sweep at modest depth to scan for the absorbing-state transition.
    local q
    for q in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50; do
        local tag=${q//./}; tag=${tag:0:3}
        run "runs/bernoulli/N24_M10000_q${q}_bern" \
            --depth 24 --prob 0.5 --q "$q" --leaves bern \
            --samples 10000 --seed "0x20240${tag}"
    done
}

bern_n26() {
    # Middle ground between bern_smoke (N=24) and bern (N=28): better finite-
    # size resolution at modest cost.
    local q
    for q in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50; do
        local tag=${q//./}; tag=${tag:0:3}
        run "runs/bernoulli/N26_M10000_q${q}_bern" \
            --depth 26 --prob 0.5 --q "$q" --leaves bern \
            --samples 10000 --seed "0x20260${tag}"
    done
}

bern_n26_upper() {
    # Upper half of the q range, where we expect the absorbing-state region.
    # q=1.00 is a sanity check (all leaves = 0 trivially gives X_N = 0).
    local q
    for q in 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00; do
        local tag=${q//./}; tag=${tag:0:3}
        run "runs/bernoulli/N26_M10000_q${q}_bern" \
            --depth 26 --prob 0.5 --q "$q" --leaves bern \
            --samples 10000 --seed "0x20260${tag}"
    done
}

bern() {
    # Production sweep at effective depth 28 to localise q_c. Refine the grid
    # once the smoke sweep tells us where the transition is.
    local q
    for q in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50; do
        local tag=${q//./}; tag=${tag:0:3}
        run "runs/bernoulli/N28_M10000_q${q}_bern" \
            --depth 28 --prob 0.5 --q "$q" --leaves bern \
            --samples 10000 --seed "0x20280${tag}"
    done
}

# --- non-critical p in {0.1, 0.2, 0.3, 0.4} ---------------------------------

noncritical() {
    # legacy noncritical sweep (M=1e4, depth+2 leaf shift). Superseded by
    # noncritical_K35 below; kept for reproducibility of any earlier reports.
    local p
    for p in 0.1 0.2 0.3 0.4; do
        local tag=${p//./}
        run "runs/noncritical/N28_M10000_p${p}_four" \
            --depth 28 --prob "$p" --leaves four --samples 10000 \
            --seed "0x30280${tag}"
    done
}

bern_K35() {
    # High-precision Bernoulli-variant sweep at p=1/2 across the full
    # q-range, at effective depth N=48 and M=1e5. Uses one mu_K=35 leaf
    # table per q value; simulator runs at depth D = 48 - 35 = 13.
    # ~5 s per q value, ~2 min total for all 21 grid points.
    local M=100000 D=13
    local q
    for q in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 \
             0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00; do
        local TABLE="tables/mu_K35_p0.5_q${q}.csv"
        if [[ ! -f "$TABLE" ]]; then
            echo "missing leaf table $TABLE" >&2; exit 1
        fi
        local tag=${q//./}
        run "runs/bernoulli/N48_M${M}_q${q}_K35" \
            --depth "$D" --prob 0.5 --leaves table --leaf-table "$TABLE" \
            --samples "$M" --seed "0x40480${tag}"
    done
}

bern_nstability() {
    # N-stability check for the conditional-law deformation reported in §4.
    # For each q in {0.25, 0.50, 0.75} we run M=1e5 at effective N=36 (D=1)
    # and N=60 (D=25), both with K=35 leaf table at the matching q.
    # Wall: N=36 trio is sub-second; N=60 trio dominates at ~16000 s each
    # (~13.3 h total). Schedule to overlap overnight.
    local M=100000
    local q effN D seed
    for q in 0.25 0.50 0.75; do
        local TABLE="tables/mu_K35_p0.5_q${q}.csv"
        if [[ ! -f "$TABLE" ]]; then
            echo "missing leaf table $TABLE" >&2; exit 1
        fi
        local tag=${q//./}
        # N=36: D=1
        run "runs/bernoulli/N36_M${M}_q${q}_K35" \
            --depth 1 --prob 0.5 --leaves table --leaf-table "$TABLE" \
            --samples "$M" --seed "0x40360${tag}"
        # N=60: D=25
        run "runs/bernoulli/N60_M${M}_q${q}_K35" \
            --depth 25 --prob 0.5 --leaves table --leaf-table "$TABLE" \
            --samples "$M" --seed "0x40600${tag}"
    done
}

noncritical_K35() {
    # High-precision noncritical sweep at effective depth N=48, M=1e5,
    # using the precomputed mu_K=35 leaf table at the corresponding p.
    # D = N - K = 48 - 35 = 13, so per-sample work ~25 * 2^13 = 2e5 ops;
    # total wall ~2 s per p value. Geometric convergence to the limit
    # P(X_inf=1) = (1-2p)/(1-p) at rate (2p)^N puts us > 5 orders of
    # magnitude inside the asymptotic regime at this depth.
    local M=100000
    local p
    for p in 0.10 0.20 0.30 0.40; do
        local TABLE="tables/mu_K35_p${p}_q0.csv"
        if [[ ! -f "$TABLE" ]]; then
            echo "missing leaf table $TABLE" >&2; exit 1
        fi
        local tag=${p//./}
        run "runs/noncritical/N48_M${M}_p${p}_K35" \
            --depth 13 --prob "$p" --leaves table --leaf-table "$TABLE" \
            --samples "$M" --seed "0x30480${tag}"
    done
}

# --- composites -------------------------------------------------------------

all_short() { smoke; bern_smoke; noncritical; }

main() {
    local cmd=${1:-smoke}
    case "$cmd" in
        smoke|medium|bern_smoke|bern_n26|bern_n26_upper|bern|bern_K35|bern_nstability|noncritical|noncritical_K35|all_short|precision_step4) "$cmd" ;;
        long)  long_run ;;
        xlong)
            echo "xlong has been disabled: superseded by precision_step4 (eff 36+ at M=1e5 via mu_K=35 leaf table)" >&2
            exit 1
            ;;
        *) sed -n '3,/^set -euo/p' "$0" | sed -n '/^# Usage/,/^#$/p' | sed 's/^# //'; exit 1 ;;
    esac
}

main "$@"
