#!/usr/bin/env bash
#
# Near-critical pilot sweep for the JPA revision (Referee 1, Comment 1.2):
# p_N = 1/2 - c/N for c in {1/2, 1, 2} and effective depths
# N in {24, 32, 40, 48, 56, 60}, M = 1e5 samples per cell, all-ones leaves (q = 0).
#
# Effective N = D + K with the usual two-table split:
#   N in {24, 32}:         K = 20 (D = 4, 12)
#   N in {40, 48, 56, 60}: K = 35 (D = 5, 13, 21, 25)
# One leaf table per (K, p_N); cells with repeated p still get distinct tables
# because K differs. Doubling cross-checks (same p_N at two depths):
#   (c=1,   N=24) and (c=2, N=48) share p = 11/24 = 0.458333...
#   (c=1/2, N=24) and (c=1, N=48) share p = 23/48 = 0.479166...
# These pairs test whether the behavior depends on c = N(1/2 - p_N) or on p alone.
#
# Seed convention: 0x70{NN}{cc}, NN = effective depth, cc in {05,10,20} for
# c in {1/2, 1, 2}.
#
# Wall estimate: the N <= 48 cells (D <= 13) run in seconds each after their
# precomputes (~5-12 min end-to-end for that block). The deep extension is the
# cost: N=56 (D=21) ~15-20 min per cell, N=60 (D=25) ~4-5 h per cell, so the
# full 18-cell grid is an overnight run. Idempotent: completed cells skip.

set -euo pipefail
cd "$(dirname "$0")"

BIN=./minplus
PRE=./precompute_leaves
test -x "$BIN" || { echo "build first: make"; exit 1; }
test -x "$PRE" || { echo "build first: make"; exit 1; }
mkdir -p runs/nearcritical logs tables

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

pre() {
    # pre <table> <p> <K>
    local table="$1" p="$2" K="$3"
    local logf="logs/precompute_$(basename "$table" .csv).log"
    if [[ -f "$table" ]]; then
        echo "[skip] ${table}"
        return 0
    fi
    echo "[pre ] ${table}"
    local t0 t1
    t0=$(date +%s)
    "$PRE" --prob "$p" --q 0 --levels "$K" --buffer-bits 22 --out "$table" >"$logf" 2>&1
    t1=$(date +%s)
    printf "[done] %s  (%ds)\n" "$table" "$((t1 - t0))"
}

M=100000

# Cells: N  ctag  p (full precision, p_N = 1/2 - c/N)  ptag (filename tag)
CELLS="
24 05 0.4791666666666667 0.4791666667
24 10 0.4583333333333333 0.4583333333
24 20 0.4166666666666667 0.4166666667
32 05 0.484375 0.484375
32 10 0.46875 0.46875
32 20 0.4375 0.4375
40 05 0.4875 0.4875
40 10 0.475 0.475
40 20 0.45 0.45
48 05 0.4895833333333333 0.4895833333
48 10 0.4791666666666667 0.4791666667
48 20 0.4583333333333333 0.4583333333
56 05 0.4910714285714286 0.4910714286
56 10 0.4821428571428571 0.4821428571
56 20 0.4642857142857143 0.4642857143
60 05 0.4916666666666667 0.4916666667
60 10 0.4833333333333333 0.4833333333
60 20 0.4666666666666667 0.4666666667
"

while read -r N ctag p ptag; do
    [[ -z "${N:-}" ]] && continue
    if (( N <= 32 )); then K=20; else K=35; fi
    D=$((N - K))
    TABLE="tables/mu_K${K}_p${ptag}_q0.csv"
    pre "$TABLE" "$p" "$K"
    run "runs/nearcritical/N${N}_M${M}_p${ptag}_c${ctag}_K${K}" \
        --depth "$D" --prob "$p" --leaves table --leaf-table "$TABLE" \
        --samples "$M" --seed "0x70${N}${ctag}"
done <<< "$CELLS"

echo "nearcritical pilot complete."
