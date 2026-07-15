#!/usr/bin/env bash
#
# One-shot precompute of every leaf table needed for the M_buffer = 22 rerun.
#
# Generates 27 tables in tables/ at K-cutoff and (p, q) values consumed by
# run_production.sh:
#   * mu_K20_p0.5_q0.csv                 -- precision_step4 shallow (eff N=24,28,32)
#   * mu_K35_p0.5_q0.csv                 -- precision_step4 deep    (eff N>=36)
#   * mu_K35_p0.5_q{0.00..1.00}.csv      -- bern_K35 sweep (21 q values)
#   * mu_K35_p{0.10,0.20,0.30,0.40}_q0.csv -- noncritical_K35 sweep (4 p values)
#
# All tables are computed with --buffer-bits 22, so j_max = 2^22 - 1 ~= 4.2e6.
# That is 4x the M=20 cap that the M=21 sensitivity test showed reduced the
# truncated mass to FFT roundoff (1.3e-15); M=22 leaves another factor of 2 of
# headroom. Each precompute is ~10-30s; total ~10 min.
#
# ONE-SHOT: refuses to overwrite existing tables. Delete tables/*.csv first if
# you really mean to regenerate.
#
# Usage:
#   ./precompute_all_M22.sh
#
# Output: tables/ populated, log to logs/precompute_M22.log.

set -euo pipefail
cd "$(dirname "$0")"

BIN=./precompute_leaves
test -x "$BIN" || { echo "build first: make precompute_leaves" >&2; exit 1; }

mkdir -p tables logs
LOGF=logs/precompute_M22.log
: > "$LOGF"

generate() {
    # generate <out> <args...>
    local out="$1"; shift
    if [[ -f "$out" ]]; then
        echo "[skip] $out (already exists)" | tee -a "$LOGF"
        return 0
    fi
    local t0 t1
    t0=$(date +%s)
    echo "[run ] $out  $*" | tee -a "$LOGF"
    "$BIN" "$@" --buffer-bits 22 --out "$out" >>"$LOGF" 2>&1
    t1=$(date +%s)
    echo "[done] $out  ($((t1 - t0))s)" | tee -a "$LOGF"
}

echo "=== K=20 critical (precision_step4 shallow) ==="
generate tables/mu_K20_p0.5_q0.csv --prob 0.5 --q 0 --levels 20

echo "=== K=35 critical (precision_step4 deep) ==="
generate tables/mu_K35_p0.5_q0.csv --prob 0.5 --q 0 --levels 35

echo "=== K=35 Bernoulli sweep (21 q values) ==="
for q in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 \
         0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00; do
    generate "tables/mu_K35_p0.5_q${q}.csv" --prob 0.5 --q "$q" --levels 35
done

echo "=== K=35 noncritical (4 p values) ==="
for p in 0.10 0.20 0.30 0.40; do
    generate "tables/mu_K35_p${p}_q0.csv" --prob "$p" --q 0 --levels 35
done

echo "=== precompute summary ==="
ls -la tables/ | tee -a "$LOGF"
echo "wrote logs/precompute_M22.log"
