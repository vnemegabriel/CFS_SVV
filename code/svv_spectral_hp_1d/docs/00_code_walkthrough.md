# 0. Code walkthrough: from the equations to the lines that compute them

`docs/01`–`04` explain the mathematics and the results. This note explains the code:
where each symbol lives, what one time step does line by line, where the Jacobian
enters, what the tests guarantee, and the C++ needed to read it. Read it first if you
have not written C++ before; read it again before changing anything in `src/`.

## 0.1 Build, test, run

Prerequisites: a C++17 compiler (g++ ≥ 9), CMake ≥ 3.16, Python 3 with numpy and
matplotlib (plots only).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release   # configure (once)
cmake --build build -j                           # compile (after every edit)
./build/test_svv_operator                        # 31 checks; must print "0 failed"
./build/burgers_svv                              # Fig. 6.27 case, writes results/
bash scripts/run_all.sh                          # everything, ~6 min
```

Run the executables from the project root: they write to `results/`. On OneDrive/WSL,
if `make` reports clock skew or does not recompile a file you edited, run
`touch src/*.cpp apps/*.cpp tests/*.cpp`.

CMake builds three kinds of targets (`CMakeLists.txt`): `polylib` (Karniadakis's
polynomial library), `svvhp` (everything in `src/`, compiled once into a static
library) and three executables (`burgers_svv`, `cns1d_svv`, `test_svv_operator`) that
link against it. A change in `src/` recompiles the library and relinks all three.

## 0.2 The four vector spaces of the code

Every `Vector` (an alias for `std::vector<double>`) in the code is one of four kinds.
Mixing them up is the most common error; the size tells them apart.

| kind | size | index | symbol | produced by |
|---|---|---|---|---|
| global modal coefficients | `A.nGlobal` ($N_{el}P + 1$, or $N_{el}P$ periodic) | global dof | $\hat u$ | `l2Projection`, time integrator |
| element modal coefficients | `nModes` $= P+1$ | $p = 0..P$ (0 and $P$ are vertices) | $\hat u_e$ | `A.gather(uhat, e)` |
| values at quadrature points | `Q` | $q = 0..Q-1$ | $u_q$ | `basis.evaluate(ue)` |
| orthonormal (Legendre) coefficients | `nModes` | degree $k = 0..P$ | $\tilde u$ | `basis.toOrthonormal(ue)` $= T \hat u_e$ |

Element and global vectors both hold modal coefficients but in different numberings:
`gather` copies global → element, `scatterAdd` adds element → global. The shared vertex
of two elements is one global entry, so `scatterAdd` *adds* the two contributions (this
is the assembly).

In `cns1d_svv` the global state stacks the three conserved variables in one `Vector`
of size `3 * nGlobal`: `[ρ̂ | m̂ | Ê]`. `solver.block(U, k)` returns a pointer to the
start of block `k`; `A.gather(ptr, e)` accepts that pointer directly.

Matrices (`Matrix`) are dense and row-major: `A(i, j)` is row `i`, column `j`.

## 0.3 Symbols and where they live

| math | code | space / scaling | defined in |
|---|---|---|---|
| $\xi \in [-1,1]$, $\xi_q$, $w_q$ | `basis.zq`, `basis.wq` | reference | `02_Basis1D` |
| $\psi_p(\xi_q)$, $\psi_p'(\xi_q)$ | `basis.B(q,p)`, `basis.dB(q,p)` | reference ($d/d\xi$) | `02_Basis1D` |
| $\tilde\psi_k(\xi_q)$ | `basis.Bortho(q,k)` | reference | `02_Basis1D` |
| $M, S, K$ | `basis.M`, `basis.S`, `basis.K` | reference (no $J$) | `02_Basis1D` |
| $T$, $T^{-1}$ | `basis.T`, `basis.Tinv` | reference | `02_Basis1D` |
| $J_e = h_e/2$, $x(\xi)$ | `mesh.J(e)`, `mesh.x(e, xi)` | — | `04_Mesh1D` |
| global mass $\sum_e J_e M$ | `A.massMatrix()` | **physical** | `05_CGAssembly` |
| global Laplacian $\sum_e K/J_e$ | `A.laplacianMatrix()` | **physical** | `05_CGAssembly` |
| $\sum_e L/J_e$ | `A.operatorMatrix(Lref)` | **physical** | `05_CGAssembly` |
| $\mathrm{dof}(e,p) = eP + p$ | `A.dof(e, p)` | — | `05_CGAssembly` |
| $\hat F_k$ | `svvKernel(P, Mcut, type)` | — | `03_SVV` |
| $\Phi = T^{-1} F T$ | `modalFilter(basis, F)` | reference | `03_SVV` |
| $L_{svv}$ (book or Kirby form) | `svvOperator(basis, F, form)`, stored as `Lref_` | reference | `03_SVV` |
| $\epsilon$ | `opt.eps` (Burgers), `svv_.eps` (CNS) | physical | `08`, `09` |
| $\epsilon_e = e\,h_e \max(\lvert u\rvert+c)/P$ | computed inside `rhs` / `applySplitSVV` when `--scaled` | physical | `09` |
| $c = \Delta t\,\epsilon/J^2$ | argument of `SVVBubbleFilter::apply` | reference mass on the left | `03`, `08`, `09` |
| $\rho(M^{-1}A)$ | `A.spectralRadius(Aglobal)` | physical | `05_CGAssembly` |

### Where the Jacobian goes

Every operator in `Basis1D` and `03_SVV` is built on the reference element. The
Jacobian $J = h/2$ enters only when an element contribution is combined with others
or with time. The rules (derived in `docs/02` §2.2):

| contribution | factor | code |
|---|---|---|
| mass $(\psi_i, u)$ | $\times J$ | `basis.M.scaled(mesh.J(e))` |
| $(\psi_i, f_x)$, $(\psi_i', G)$ | none (a $J$ and a $1/J$ cancel) | `b.S * b.project(F)`, `b.innerProductDeriv(G)` |
| gradient at a point, $u_x$ | $\times 1/J$ | `scaled(b.evaluateDeriv(uh), 1.0 / J)` in `09` |
| $(\psi_i', u_x)$, $(\psi_i', F * u_x)$ | $\times 1/J$ | `axpy(-eps / J, Lref_ * ue, re)` |
| any element-local solve with the *reference* mass | divide the whole equation by $J$ | $c = \Delta t\,\epsilon/J^2$ |

When you add a term, write it first as a weak form in $x$, change variables to $\xi$,
and count the powers of $J$ before writing code. A missing factor of $J$ passes
every operator-level test (symmetry, null space, $F = I \Rightarrow K$) and only shows up
when the element is compared with its physical system. That is how the split
application ran with $\epsilon J$ instead of $\epsilon$ (`docs/01` §1.4).

## 0.4 One time step of `burgers_svv`, line by line

The loop in `apps/burgers_svv.cpp`:

```cpp
RHSFunction R = [&](const Vector &u, double t) { return solver.rhs(u, t); };
for (int n = 0; n < nsteps; ++n)
{
    rk4Step(uhat, t, dt, R);          // 1. explicit convection (+ Galerkin SVV)
    solver.applySplitSVV(uhat, dt);   // 2. implicit SVV on bubbles (split only)
    t += dt;
}
```

**1. `rk4Step`** (`06_TimeIntegration.h`) calls `R` four times. Each call is
`BurgersSolver::rhs(uhat, t)` (`08_Burgers.cpp`):

| line (abridged) | what it computes | space |
|---|---|---|
| `Vector r(A_.nGlobal, 0.0);` | global residual, zeroed | global modal |
| `ue = A_.gather(uhat, e);` | coefficients of element `e` | element modal |
| `uq = b.evaluate(ue);` | $u(\xi_q) = \sum_p \hat u_p \psi_p(\xi_q)$ | quadrature |
| `fq[q] = 0.5 * uq[q] * uq[q];` | $f(u)$ at the points (aliased if `Q` small) | quadrature |
| `b.project(fq)` | $\hat f = M^{-1} B^\top W f_q$ | element modal |
| `re = scaled(b.S * ..., -1.0);` | $-(\psi_i, f_x) = -S\hat f$ | element modal |
| `axpy(-eps/J, Lref_ * ue, re);` | $-(\epsilon/J) L \hat u_e$ (Galerkin only) | element modal |
| `A_.scatterAdd(re, e, r);` | assemble | global modal |
| `A_.zeroDirichlet(r);` | $\dot u = 0$ at Dirichlet dofs | global modal |
| `return massLU_.solve(r);` | $\dot{\hat u} = M^{-1} r$ (LU factored once in the constructor) | global modal |

**2. `applySplitSVV`** (split only): for each element, `gather`, call
`bubbleFilter_.apply(ue, dt * eps / (J * J))`, and write back only the bubble
entries `p = 1..P-1`. The vertices are never written, so continuity is untouched.
`SVVBubbleFilter` factorises $M_{bb} + c L_{bb}$ once and refactorises only when
`c` changes.

### Differences in `cns1d_svv`

`CompressibleNS1D::rhs` (`09_CompressibleNS.cpp`) runs the same loop for three
variables at once:

1. `gather` $\hat\rho_e, \hat m_e, \hat E_e$ and evaluate them and their $x$-gradients
   at the quadrature points.
2. `pointState` turns $(\rho, m, E)$ and their gradients into $u, p, T, u_x, T_x$ by the
   chain rule, and **throws** if $\rho \le 0$ or $p \le 0$. The app catches the
   exception, prints `ABORT`, and still writes the last state.
3. Fluxes $F$ and $G$ at the points; $-S\hat f$ for convection, $-B'^\top W G_q$ for
   viscosity, $-(\epsilon_e/J) L \hat U_e$ for Galerkin SVV.
4. One `scatterAdd` per variable into its block; three mass solves at the end.

The time step is recomputed every 10 steps (`stableTimeStep`), and `applySplitSVV`
filters each variable (density optional, `--nosvv-rho`).

## 0.5 The C++ you need for this repository

Only the constructs that appear in the code, with the place to see each one.

| construct | example in the repo | what to know |
|---|---|---|
| header / source split | `03_SVV.h` declares, `03_SVV.cpp` defines | `#pragma once` prevents double inclusion; the `.h` is the interface, read it first |
| type alias | `using Vector = std::vector<double>;` (`01_DenseMatrix.h`) | `Vector` is a plain `std::vector`: `.size()`, `[i]`, `.data()` |
| `const` reference argument | `Vector rhs(const Vector &uhat, ...)` | no copy, read-only. Without `&` the vector is copied on every call |
| `const` member function | `Vector rhs(...) const;` | does not modify the object. `applySplitSVV` is not `const` because the filter caches its LU |
| operator overloading | `b.S * b.project(fq)`, `A(i, j)` | `*` is matrix product, `()` is element access |
| `struct` with defaults | `BurgersOptions`, `GasProperties` | fields have default values; the app overwrites them from the command line |
| `enum class` | `SVVApplication::Split` | a named choice; `parseSVVApplication("split")` converts the CLI string |
| constructor initializer list | `BurgersSolver::BurgersSolver(...) : A_(A), opt_(opt), kernel_(...)` | **members are initialised in declaration order, not list order.** `bubbleFilter_` needs `Lref_`, so it is declared last in `08_Burgers.h` |
| lambda | `[&](const Vector &u, double t) { return solver.rhs(u, t); }` | an inline function; `[&]` captures local variables by reference |
| `std::function` | `using RHSFunction = std::function<Vector(const Vector &, double)>;` | a variable that holds any callable with that signature; lets `rk4Step` work for Burgers and CNS |
| `std::unique_ptr` | `std::unique_ptr<LUSolver> lu_;` in `SVVBubbleFilter` | owns an object created later (`std::make_unique`); freed automatically |
| exceptions | `throw std::runtime_error(...)` in `pointState`, `try/catch` in `cns1d_svv.cpp` | an error jumps out of all nested calls to the nearest `catch` |
| raw pointer into a vector | `double *block(Vector &U, int k)` returns `U.data() + k * nGlobal` | pointer arithmetic: points at element `k * nGlobal`; no bounds check |
| reference member | `const Basis1D &basis_;` in `SVVBubbleFilter` | the object must outlive the filter; here `basis` lives in `main` |

Numerical traps specific to C++:

* **Integer division.** `1/16` is `0`. Write `1.0 / 16.0`. `P / 2` with `int P = 15` is `7`.
* **Command-line numbers.** `--eps 1/64` is parsed by `std::stod` as `1`. Pass
  decimals: `--eps 0.015625`.
* **No bounds checking.** `A(i, j)` and `v[i]` do not check the index. An out-of-range
  index gives garbage or a crash far from the cause. For debugging, build with
  `-DCMAKE_BUILD_TYPE=Debug` and add `-fsanitize=address,undefined` to the compile
  options.
* **Floating-point comparison.** Tests use tolerances (`< 1e-9`), never `==`, except
  where the value is exactly assigned (the filter's vertex coefficients).

## 0.6 What the tests guarantee, and what they do not

`tests/test_svv_operator.cpp` (31 checks) covers:

* **Basis**: $T^\top T = M$, $T^{-1} = M^{-1}T^\top$, $S^\top M^{-1} S = K$, symmetry,
  vertex modes are the only non-zero modes at $\xi = \pm 1$, exact differentiation.
* **Kernel**: $\hat F_k = 0$ for $k \le M_{SVV}$, $\hat F_P = 1$, monotone.
* **Operator**: both forms symmetric PSD; $F = I$ gives $K$; null spaces (degree
  $\le M_{SVV}+1$ and $\le M_{SVV}$); $L \le K$ in energy; $M\Phi$ symmetric; vertex
  rows and columns of $L$ zero.
* **Split filter**: leaves resolved polynomials unchanged, damps the top mode, never
  touches vertices, dissipative in $L^2$, and **one split step equals implicit Euler
  of the physical element system** (the $J$ scaling).
* **Assembly**: dof counts, shared vertices, periodic wrap, exact representation of
  $u = 1$ and exact projection of a cubic, range of $\rho(M^{-1}K)$.

Not covered by unit tests (covered only by the application runs in `docs/03`–`04`):

* the convective residuals (checked by the pre-shock Burgers error, `docs/03` §3.1);
* the Galerkin SVV residual scaling $(\epsilon/J)L$ against the physical system;
* the compressible RHS, `pointState` and the viscous fluxes (checked by the entropy
  wave and the viscous Sod tube);
* the time integrators and the time-step estimates.

A change to any of these needs a run of the corresponding case, not only the tests.

## 0.7 Command-line reference

`burgers_svv`:

| option | default | meaning |
|---|---|---|
| `--nel N` | 5 | elements on $[-1,1]$ |
| `--P P` | 15 | polynomial order (P+1 modes) |
| `--Q Q` | P+1 | GLL quadrature points for the non-linear term |
| `--T T` | 0.5 | final time |
| `--cfl c` | 0.5 | advective CFL, $\Delta t = c\,h_{min}/(1.2\max\lvert u_0\rvert)$ |
| `--dt dt` | — | override the time step |
| `--ic sine\|step` | sine | initial condition |
| `--conv cons\|nc\|skew` | cons | convective form |
| `--nosvv` | off | disable SVV |
| `--eps e` | 0.0625 | physical SVV amplitude |
| `--Mcut M` | 8 | kernel cut-off $M_{SVV}$ |
| `--kernel exp\|step` | exp | eq. 6.5.13 or 6.5.12 |
| `--form book\|kirby` | book | $L_{svv}$ or $\Phi^\top K \Phi$ |
| `--svv-apply split\|galerkin` | split | how the operator is applied |
| `--nprint n` | 500 | print every n steps |
| `--out tag` | automatic | writes `results/tag.csv`, `results/tag_spectrum.csv` |

`cns1d_svv`:

| option | default | meaning |
|---|---|---|
| `--case entropywave\|sod\|shuosher` | sod | test case (sets domain, BCs, $T$) |
| `--nel N`, `--P P` | 40, 8 | mesh and order |
| `--Q Q` | P+2 | quadrature points |
| `--T T` | per case | final time |
| `--cfl c` | 0.4 | advective CFL |
| `--mu`, `--Pr`, `--gamma` | 0, 0.72, 1.4 | gas; `--mu 0` is Euler |
| `--ic-smooth δ` | 0 | tanh width of initial jumps |
| `--svv` | off | enable SVV |
| `--eps e` | 0.05 | constant $\epsilon$, or $e$ in $\epsilon_e = e\,h_e\max(\lvert u\rvert+c)/P$ with `--scaled` |
| `--scaled` | off | use the $\mathcal{O}(h/P)$ amplitude law |
| `--Mcut M` | P/2 | kernel cut-off |
| `--kernel`, `--form`, `--svv-apply` | exp, book, split | as in Burgers |
| `--nosvv-rho` | off | do not filter density |
| `--nprint n`, `--out tag` | 200, automatic | output control |

## 0.8 How to extend

**A new compressible case.** In `apps/cns1d_svv.cpp`, add an `else if (cas == "name")`
branch that sets `xa, xb, T, bc` and the lambdas `rho0, u0, p0` (and `rhoExact` if
known). Nothing in `src/` changes.

**A new SVV kernel.** In `03_SVV.h` add a value to `SVVKernelType`; in `03_SVV.cpp`
add its string to `parseKernelType` and a `case` in `svvKernel`. Add a test that the
kernel is zero for $k \le M_{SVV}$ and one at $k = P$, and rerun the Burgers sweep.

**A new term in the RHS.** Write the weak form in $x$, count the powers of $J$
(§0.3), add the element contribution before `scatterAdd`, and add a test that compares
the element contribution with a hand-computed physical value on an element with
$J \ne 1$. Tests on $J = 1$ cannot see a missing Jacobian.

**A new test.** Use the `CHECK(condition, "message")` macro in
`tests/test_svv_operator.cpp`. A test is useful only if it fails on the bug it
targets: break the code on purpose once and check that the test fails.

## 0.9 Debugging checklist

1. `./build/test_svv_operator` first. If a check fails, the problem is in `02`–`05`.
2. Read the start-up print: kernel values, $h_{min}$, $\Delta t$, spectral radius. A
   wrong option shows up there.
3. Reduce the problem: `--nel 1`, small `--P`, short `--T`, `--nprint 1`.
4. Time-step dependence: halve `--dt`. The Galerkin results must not change; split
   results change at first order (`docs/03` §3.3).
5. Look at the spectra (`*_spectrum.csv`): which modes grow, in which element.
6. An `ABORT` in `cns1d_svv` reports $x$, $t$, $\rho$, $p$; the CSV still holds the last
   valid state for plotting.
