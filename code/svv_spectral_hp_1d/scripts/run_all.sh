#!/usr/bin/env bash
# Build, test, run every case, and make the figures. Run from WSL:
#   cd /mnt/c/Users/valne/OneDrive/TESIS/code/svv_spectral_hp_1d && bash scripts/run_all.sh
set -uo pipefail
cd "$(dirname "$0")/.."

touch src/*.cpp apps/*.cpp tests/*.cpp          # beat OneDrive/WSL clock skew
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release > /dev/null
cmake --build build -j || exit 1
mkdir -p results

echo "=== unit tests ==="
./build/test_svv_operator || exit 1

echo "=== Burgers, K&S Fig. 6.27 (5 elements x 16 modes, M_SVV=8, eps=1/16) ==="
./build/burgers_svv --nosvv                      # (a)
./build/burgers_svv                              # (b) split application, book parameters (eps = 1/16)
./build/burgers_svv --svv-apply galerkin         # (b') weak-form residual: interface spikes
python3 scripts/plot_burgers.py
python3 scripts/plot_compare.py results/burgers_galerkin_vs_split.png burgers_svv burgers_svv_galerkin
./build/burgers_svv --eps 0.0125 --out burgers_svv_e80                          # eps = 1/80
./build/burgers_svv --eps 0.0125 --svv-apply galerkin --out burgers_svv_galerkin_e80
python3 scripts/plot_compare.py results/burgers_galerkin_vs_split_e80.png burgers_svv_e80 burgers_svv_galerkin_e80

echo "=== Burgers parameter sweep ==="
bash scripts/sweep_burgers.sh

echo "=== Compressible: entropy wave (spectral accuracy with SVV on) ==="
for P in 4 6 8 10 12; do
  printf 'P=%2d  no SVV: ' $P; ./build/cns1d_svv --case entropywave --nel 8 --P $P --out ew_nosvv | grep 'L2 error'
  printf 'P=%2d     SVV: ' $P; ./build/cns1d_svv --case entropywave --nel 8 --P $P --svv --eps 0.05 --scaled --out ew_svv | grep 'L2 error'
done
python3 scripts/plot_cns.py --case entropywave results/ew_nosvv.csv results/ew_svv.csv --out results/cns_entropywave.png

echo "=== Compressible: viscous Sod shock tube (mu = 2e-3 and 5e-4) ==="
./build/cns1d_svv --case sod --nel 40 --P 8 --mu 2e-3 --out sod_ns_nosvv | tail -3
./build/cns1d_svv --case sod --nel 40 --P 8 --mu 2e-3 --svv --eps 0.1 --scaled --out sod_ns_svv | tail -3
python3 scripts/plot_cns.py --case sod results/sod_ns_nosvv.csv results/sod_ns_svv.csv --out results/cns_sod_ns.png
./build/cns1d_svv --case sod --nel 40 --P 8 --mu 5e-4 --out sod_ns2_nosvv | tail -3
./build/cns1d_svv --case sod --nel 40 --P 8 --mu 5e-4 --svv --eps 0.1 --scaled --out sod_ns2_svv | tail -3
python3 scripts/plot_cns.py --case sod results/sod_ns2_nosvv.csv results/sod_ns2_svv.csv --out results/cns_sod_ns2.png

echo "=== Compressible: inviscid Sod, M_SVV = P/2 (expected to FAIL) ==="
./build/cns1d_svv --case sod --nel 40 --P 8 --ic-smooth 0.02 --svv --eps 0.3 --scaled --out sod_euler_svv | grep -E 'ABORT|final'
echo "=== Compressible: inviscid Sod with M_SVV = 1 (SVV on almost every mode: survives) ==="
./build/cns1d_svv --case sod --nel 40 --P 8 --Mcut 1 --ic-smooth 0.02 --svv --eps 0.3 --scaled --out sod_euler_svv_m1 | grep -E 'ABORT|final'
python3 scripts/plot_cns.py --case sod results/sod_euler_svv_m1.csv --out results/cns_sod_euler_m1.png

echo "=== Compressible: Shu-Osher, viscous (mu = 2e-2) ==="
./build/cns1d_svv --case shuosher --nel 100 --P 8 --mu 2e-2 --ic-smooth 0.1 --out so_nosvv | tail -2
./build/cns1d_svv --case shuosher --nel 100 --P 8 --mu 2e-2 --ic-smooth 0.1 --svv --eps 0.1 --scaled --out so_svv | tail -2
python3 scripts/plot_cns.py --case shuosher results/so_nosvv.csv results/so_svv.csv --out results/cns_shuosher.png
echo "done: figures are in results/*.png"
