#!/usr/bin/env bash
# Parameter sweep for the Burgers test. Prints TV and error for each variant.
#   bash scripts/sweep_burgers.sh
cd "$(dirname "$0")/.."
run() { tag=$1; shift; printf '%-42s' "$*"; ./build/burgers_svv "$@" --out "$tag" | grep -E 'BLOW|total variation' | sed 's/(exact.*//; s/ *total variation of u_h at T: /TV=/' | tr '\n' ' '; python3 scripts/burgers_error.py "results/$tag.csv" | sed 's/results.[a-z0-9_]*.csv: //'; }
echo "--- no SVV"
run sw --nosvv --Q 16
run sw --nosvv --Q 24
run sw --nosvv --Q 24 --conv nc
run sw --nosvv --Q 16 --conv skew
echo "--- SVV, split application (element-wise implicit bubble filter)"
run sw --Q 16
run sw --Q 17
run sw --Q 24
run sw --Q 24 --conv nc
run sw --Q 16 --conv skew
run sw --Q 16 --form kirby
run sw --Q 16 --Mcut 4
run sw --Q 16 --Mcut 12
run sw --Q 16 --eps 0.25
run sw --Q 16 --eps 0.015625
run sw --Q 16 --eps 0.0125
run sw --Q 16 --eps 0.0125 --Mcut 4
run sw --Q 16 --eps 0.00390625
run sw --Q 16 --kernel step
echo "--- SVV, Galerkin application (weak-form residual through the global mass matrix)"
run sw --Q 16 --svv-apply galerkin
run sw --Q 24 --svv-apply galerkin
run sw --Q 16 --svv-apply galerkin --eps 0.25
run sw --Q 16 --svv-apply galerkin --eps 0.0125
run sw --Q 16 --svv-apply galerkin --eps 0.00390625
