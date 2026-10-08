# svv_spectral_hp_1d — Spectral Vanishing Viscosity in a 1D C0 spectral/hp code

A from-scratch C++17 implementation of the SVV test of Karniadakis & Sherwin,
*Spectral/hp Element Methods for CFD*, §6.5.2 (Fig. 6.27, inviscid Burgers with
continuous Galerkin), extended to the 1D compressible Euler/Navier–Stokes equations.
No dependencies beyond a C++ compiler and CMake (plus Python/matplotlib for plots).
It is written to be read: files are numbered in reading order, every header starts
with the derivation it implements, and the book's identities are unit-tested.

## Build and run (WSL Ubuntu)

```bash
cd /mnt/c/Users/valne/OneDrive/TESIS/code/svv_spectral_hp_1d
bash scripts/run_all.sh          # build, unit tests, every case, every figure
```

or step by step:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
./build/test_svv_operator                       # 31 checks of the book's identities and of the scaling
./build/burgers_svv --nosvv                     # Fig. 6.27(a)
./build/burgers_svv                             # Fig. 6.27(b)
python3 scripts/plot_burgers.py                 # -> results/fig_6_27_burgers.png
./build/cns1d_svv --case sod --nel 40 --P 8 --mu 2e-3 --svv --eps 0.1 --scaled
./build/burgers_svv --help ; ./build/cns1d_svv --help
```

Run the executables from the project root: they write to `results/`. If WSL reports
"clock skew" (OneDrive time stamps), `touch src/*.cpp apps/*.cpp` before building;
`run_all.sh` does it for you.

## Layout and reading order

New to the code (or to C++)? Start with `docs/00_code_walkthrough.md`.


```
extern/polylib/      Karniadakis' polynomial library (same code as in Nektar++)
src/01_DenseMatrix   dense matrix + LU (so nothing is hidden in a library)
src/02_Basis1D       modified C0 basis, orthonormal Legendre, M S K T matrices
src/03_SVV           kernel (6.5.13), operator (p. 351), modal filter, bubble filter
src/04_Mesh1D        mesh, Jacobian
src/05_CGAssembly    local->global map, assembly, BCs, L2 projection, spectral radius
src/06_TimeIntegration  RK4 / SSP-RK3
src/07_Output        sampling, modal spectra, CSV, CLI args
src/08_Burgers       Burgers RHS (three convective forms) + SVV
src/09_CompressibleNS   1D compressible NS RHS + SVV
apps/burgers_svv.cpp    the Fig. 6.27 test
apps/cns1d_svv.cpp      entropy wave / Sod / Shu-Osher
tests/test_svv_operator.cpp
scripts/   plot_burgers.py, plot_cns.py, exact_riemann.py (Toro), sweeps
docs/      00 code walkthrough  01 theory  02 discretisation  03 Burgers results  04 compressible results
```

Suggested path: `docs/00` → `docs/01` → `src/02`, `src/03` and the tests → `docs/02` → `src/05`,
`src/08`, `apps/burgers_svv.cpp` → `docs/03` and play with the sweep →
`src/09`, `apps/cns1d_svv.cpp`, `docs/04`.

## What was found (details and tables in `docs/03` and `docs/04`)

* **The operator of p. 351 is verified**: `TᵀT = M`, `T⁻¹ = M⁻¹Tᵀ`, `SᵀM⁻¹S = K`,
  `L_svv` symmetric PSD, `F = I` gives the Laplacian, and `L_svv` annihilates every
  polynomial of degree ≤ `M_SVV+1` (no viscosity on resolved modes).
* **Scaling bug fixed (September 2026)**: the split application used `c = Δt ε/J`
  instead of `Δt ε/J²`, i.e. it applied `εJ` instead of `ε`. Every result produced
  before the fix ran with 5× (Burgers) to 80× (Sod) less SVV than stated. A unit test
  now compares one split step with the physical element system.
* **Fig. 6.27(b) is not reproduced with the book's ε = 1/16** taken as a physical
  viscosity: both applications raise TV above the no-SVV run (split 10.9, Galerkin
  15.8, no SVV 7.4; exact 4). The figure's look appears at ε ≈ 1/80 (TV 6.8). The
  book's ε convention remains open (`docs/03` §3.4).
* **The vertex rows and columns of the SVV operator are zero.** Through the consistent
  mass matrix this makes the Galerkin application worse than the split one at equal
  ε, but the split application also shows interface spikes at the book's amplitude.
  The split application is a first-order splitting: its result depends on Δt
  (`docs/03` §3.3).
* **SVV bounds aliasing-unstable convective forms** (TV 55 → 13 for the
  non-conservative form) without making them better than the stable conservative form.
* **Compressible**: spectral convergence is preserved with SVV on, at a cost of one to
  two orders of magnitude in error for P ≤ 8 (entropy wave); on the viscous Sod tube
  SVV brings TV(ρ) from 1.80 to the exact 0.875; the inviscid Sod tube fails with the
  default `M_SVV = P/2` but survives with `M_SVV ≤ 2` (`cns_sod_euler_m1.png`);
  Shu–Osher needs μ = 2e-2 with or without SVV.

## Figures

| file | content |
|---|---|
| `results/fig_6_27_burgers.png` | Fig. 6.27 reproduction (a) no SVV (b) SVV |
| `results/burgers_spectra.png` | per-element Legendre spectra with/without SVV |
| `results/burgers_galerkin_vs_split.png` | the two ways of applying the operator, ε = 1/16 |
| `results/burgers_galerkin_vs_split_e80.png` | same, ε = 1/80 |
| `results/cns_entropywave.png` | smooth verification case |
| `results/cns_sod_ns.png`, `cns_sod_ns2.png` | viscous Sod, μ = 2e-3 and 5e-4, exact overlay |
| `results/cns_sod_euler_m1.png` | inviscid Sod, M_SVV = 1 |
| `results/cns_shuosher.png` | viscous Shu–Osher, μ = 2e-2 |

## Command-line reference

`burgers_svv`: `--nel --P --Q --T --cfl --dt --nosvv --eps --Mcut --kernel exp|step
--form book|kirby --conv cons|nc|skew --svv-apply split|galerkin --ic sine|step --out`

`cns1d_svv`: `--case entropywave|sod|shuosher --nel --P --Q --T --cfl --mu --Pr --gamma
--ic-smooth --svv --eps [--scaled] --Mcut --kernel --form --svv-apply --nosvv-rho --out`
