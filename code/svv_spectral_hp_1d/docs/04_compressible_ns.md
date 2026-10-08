# 4. Extension to 1D compressible Navier–Stokes

`src/09_CompressibleNS.*`, `apps/cns1d_svv.cpp`. This is the 1D prototype of the
thesis plan's Phase 0 ("SVV on the conserved variables of the compressible solver").

> **Correction (September 2026).** Earlier versions of this note were produced with a
> split application that applied $\epsilon_e J_e$ instead of $\epsilon_e$ (see `docs/01`
> §1.4): 8 times too little SVV on the entropy wave, 80 times on Sod, 20 times on
> Shu–Osher. All numbers below use the corrected code. The main change is in §4.3:
> the inviscid Sod tube is no longer a uniform failure.

## 4.1 Equations and weak form

$$
U = (\rho,\ m = \rho u,\ E), \qquad
p = (\gamma - 1)\left(E - \frac{m^2}{2\rho}\right), \qquad
T = \frac{p}{\rho R}
$$

$$
U_t + F(U)_x = G(U, U_x)_x + \epsilon\, \partial_x\!\left[F * U_x\right]
$$

$$
F = \begin{pmatrix} m \\ m u + p \\ u(E + p) \end{pmatrix}, \qquad
G = \begin{pmatrix} 0 \\ \tau \\ u\tau - q \end{pmatrix}, \qquad
\tau = \tfrac43\, \mu\, u_x, \qquad
q = -\kappa\, T_x, \qquad
\kappa = \frac{\mu\, c_p}{\mathrm{Pr}}
$$

Each conserved variable gets its own copy of the Burgers machinery:

$$
(v, U_t) = -(v, F_x) - (v_x, G) - \epsilon\, (v_x, F * U_x)
$$

$$
\text{per element:} \qquad -S \hat f \;-\; B'^\top W G_q \;-\; \frac{\epsilon}{J}\, L_{svv}\, \hat U
$$

The same mass matrix (with Dirichlet rows) serves the three variables. The default
quadrature is $Q = P + 2$ (one point of over-integration), unlike Burgers ($Q = P+1$).

## 4.2 Choices that are specific to the compressible case

* **SVV on which variables?** Default: all three conserved variables (`--nosvv-rho`
  switches density off). Entropy-variable SVV (thesis Phase 4) would go in
  `applySplitSVV` as a change of variables before/after the filter.
* **Amplitude.** `--eps e` is a constant; with `--scaled`,
  $\epsilon_e = e\, h_e \max(|u| + c)/P$, an $\mathcal{O}(h/P)$ viscosity in the spirit of Kirby &
  Karniadakis 2002 and Nektar++'s `SVVDiffCoeff`. Because the element operator
  carries $1/J_e = 2/h_e$, the weak-form coefficient $\epsilon_e / J_e = 2\, e \max(|u| + c)/P$ is independent of
  the element size: the coefficient is a pure "fraction of a grid-scale viscosity".
  The damping *rate* is not: dividing by the element mass $J_e M$ gives
  $\epsilon_e/J_e^2 \sim 1/h$, and the split step uses $c = \Delta t\, \epsilon_e / J_e^2$.
* **Application.** `--svv-apply split` (default) or `galerkin`, as in Burgers.
* **Positivity.** The code throws as soon as $\rho \le 0$ or $p \le 0$ at a quadrature
  point, reporting where and when. SVV has no positivity property: this is the
  main limitation you will hit.
* **Time step.** Recomputed every 10 steps from $\max(|u| + c)$, the physical
  diffusivity $\max\!\left(\dfrac{4\mu}{3\rho}, \dfrac{\kappa}{\rho c_v}\right)$ times $\rho(M^{-1}K)$, and (Galerkin application only)
  the SVV spectral radius.

## 4.3 Cases and results

### Entropy wave (verification, periodic, exact solution)

$\rho = 1 + 0.2 \sin\!\left(\pi(x - t)\right)$, $u = p = 1$, 8 elements ($J = 1/8$), $T = 2$. $L^2$ error of $\rho$:

| $P$ | no SVV | SVV (`--eps 0.05 --scaled`) | ratio |
|---|---|---|---|
| 4 | $2.5 \times 10^{-6}$ | $5.3 \times 10^{-5}$ | 21 |
| 6 | $1.7 \times 10^{-9}$ | $9.7 \times 10^{-8}$ | 57 |
| 8 | $2.2 \times 10^{-11}$ | $5.9 \times 10^{-10}$ | 27 |
| 10 | $4.3 \times 10^{-12}$ | $1.6 \times 10^{-11}$ | 3.6 |
| 12 | $9.3 \times 10^{-13}$ | $9.3 \times 10^{-13}$ | 1.0 |

Spectral convergence is preserved with SVV on, but the cost at low and moderate $P$
is one to two orders of magnitude, not "a factor of a few". The kernel starts at
$M_{SVV} = P/2$, so for $P \le 8$ it acts on modes that carry the resolved wave. From
$P = 12$ on, the difference is nil. Reducing $\Delta t$ (`--cfl 0.1`) changes these
errors by less than 2 %, so the split error is negligible here. This is the "level-4
test" of the thesis testing ladder; the ratio column is the number to watch when the
amplitude law changes.

### Viscous Sod shock tube ($[0,1]$, $T = 0.2$, 40 elements, $P = 8$, `--eps 0.1 --scaled`)

| $\mu$ | $\mathrm{TV}(\rho)$ no SVV | $\mathrm{TV}(\rho)$ SVV | exact $\mathrm{TV}(\rho)$ |
|---|---|---|---|
| $2 \times 10^{-3}$ | 1.80 | 0.875 | 0.875 |
| $5 \times 10^{-4}$ | 1.70 | 0.875 | 0.875 |

Figures `cns_sod_ns.png`, `cns_sod_ns2.png`: without SVV the head of the rarefaction
carries a train of oscillations (Gibbs from the initial jump, transported and never
damped because the physical viscosity is small there); SVV removes them while the
contact and the (resolved, viscous) shock are unchanged. The density TV equals the
exact value to three digits. Exact Riemann solution is overlaid
(`scripts/exact_riemann.py`, Toro's iterative solver).

### Inviscid Sod

$\mu = 0$, 40 elements, `--scaled`, smoothed initial jump (`--ic-smooth 0.02`; without
smoothing, $P = 8$, $M_{SVV} = 1$, $e = 0.3$ fails at the first step from the $t = 0$ Gibbs oscillations):

| $P$ | $M_{SVV}$ | $e = 0.01$ | $e = 0.3$ | $e = 3$ |
|---|---|---|---|---|
| 4 | 1 | fails | $\mathrm{TV}(\rho) = 1.16$ | $\mathrm{TV}(\rho) = 1.09$ |
| 4 | 4 | fails | fails | fails |
| 8 | 1 | fails | $\mathrm{TV}(\rho) = 1.15$ | $\mathrm{TV}(\rho) = 1.38$ |
| 8 | 2 | | $T$ reached | $T$ reached |
| 8 | 3 | | fails ($t \approx 0.07$) | fails ($t \approx 0.18$) |
| 8 | 4 ($= P/2$) | fails | fails ($t \approx 0.07$) | fails ($t \approx 0.07$) |

"Fails" means negative pressure at the shock. With the default cut-off $M_{SVV} = P/2$
SVV does not survive the inviscid shock. With $M_{SVV} \le 2$ and enough amplitude
it does: `cns_sod_euler_m1.png` ($P = 8$, $M_{SVV} = 1$, $e = 0.3$) shows the correct
shock position, a smeared contact and rarefaction, and an undershoot of $u$ at the
shock (about $-0.45$). This is consistent with the argument that a shock needs an
$\mathcal{O}(h/P)$ viscosity on almost all modes: with $M_{SVV} = 1$ SVV *is* such a
viscosity, and it is no longer "vanishing". The division of labour in the thesis plan
stands (SVV for under-resolution, sensor-gated artificial viscosity for shocks), but
the negative result is narrower than earlier versions of this note stated.

### Shu–Osher ($[-5,5]$, $T = 1.8$, 100 elements, $P = 8$, `--ic-smooth 0.1`)

The Mach-3 shock is only survivable here as a *resolved viscous* shock:

| $\mu$ | no SVV | SVV (`--eps 0.1 --scaled`) |
|---|---|---|
| $2 \times 10^{-3}$ | aborts $t \approx 0.07$ | aborts $t \approx 0.08$ |
| $1 \times 10^{-2}$ | aborts $t \approx 0.14$ | aborts $t \approx 0.14$ |
| $2 \times 10^{-2}$ | runs, $\mathrm{TV}(\rho) = 9.37$ | runs, $\mathrm{TV}(\rho) = 7.38$ |

At $\mu = 2 \times 10^{-2}$ the shock is thick enough to resolve and SVV reduces the spurious
variation by about 21 % while leaving the entropy waves behind the shock intact
(`cns_shuosher.png`). The case is included as a stress test and as the starting
point for the sensor-gated SVV of the thesis.

## 4.4 What to carry to Nektar++

* The operator is exactly `StdExpansion::v_SVVLaplacianFilter`'s (projection to the
  orthogonal basis, diagonal scaling, back-projection) applied to a gradient.
* The insertion point in `CompressibleFlowSolver` is an `ArtificialDiffusion`
  plugin: replace the pointwise $\mu(x)\, \nabla U$ by the per-element modal operation on
  `qfield` (see `SVV_nektar_implementation_plan.md`, §1.3). In DG the interface
  question of `docs/01` §1.4 does not arise in the same way (no shared vertex dofs;
  coupling goes through numerical fluxes), but the explicit-time-step consequence
  does: the spectral radius of the SVV operator must enter `GetElmtTimeStep`.
* **Check the Jacobian bookkeeping first.** The bug corrected here (one factor of
  $J$ missing when an operator built on the reference element is combined with a
  reference mass matrix) is the kind that survives every operator-level test. The
  Nektar++ implementation needs a test of one SVV step against the physical element
  system, like `split SVV step = implicit Euler of the physical element system` in
  `tests/test_svv_operator.cpp`.
* The `--scaled` amplitude is the $\texttt{SVVDiffCoeff} \cdot \dfrac{h}{P}\, (|u| + c)$ law; measure its
  effect on the isentropic vortex convergence plot exactly as done here for the
  entropy wave.
