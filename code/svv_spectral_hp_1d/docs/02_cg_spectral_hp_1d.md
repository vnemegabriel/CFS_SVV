# 2. The 1D continuous-Galerkin spectral/hp code, file by file

Read the sources in numerical order. Each header starts with a derivation block;
this note is the map.

| File | What it is | K&S reference |
|---|---|---|
| `01_DenseMatrix` | tiny dense matrix / LU. No dependencies. | — |
| `02_Basis1D` | modified $C^0$ basis, orthonormal Legendre, quadrature, $M, S, K, T$ | Ch. 2, §2.3 |
| `03_SVV` | kernel, $L_{svv}$, modal filter, bubble filter | §6.5.2 |
| `04_Mesh1D` | vertices, Jacobian $J = h/2$, affine map | §4.1 |
| `05_CGAssembly` | local↔global map, assembly, BCs, projection, spectral radius | Ch. 4 |
| `06_TimeIntegration` | RK4 / SSP-RK3 | — |
| `07_Output` | sampling, spectra, CSV, argument parsing | — |
| `08_Burgers` | RHS of Burgers + SVV | §6.5.2 |
| `09_CompressibleNS` | RHS of 1D compressible NS + SVV | — |

## 2.1 Element and quadrature

Everything is done on the reference interval $\xi \in [-1,1]$. `polylib` (Karniadakis's
own library, identical to the one inside Nektar++) provides

* `zwglj`: Gauss–Lobatto–Legendre points/weights (exact to degree $2Q-3$);
* `jacobfd` / `jacobd`: Jacobi polynomials $P^{\alpha,\beta}_n$ and their derivatives.

`Basis1D` keeps two quadratures:

* an internal exact one ($Q_e = P+2$) used once to build the reference matrices;
* the *solution* quadrature ($Q \ge P+1$, user-chosen) used for the non-linear terms.

$Q = P+1$ is the classic choice (the GLL points are also the interpolation nodes of
a nodal basis, so this is "collocation") and the default of `burgers_svv`;
`cns1d_svv` defaults to $Q = P+2$. $Q > P+1$ is over-integration. The Burgers
sweep shows the two choices behave quite differently under aliasing.

Matrices, all $(P+1)\times(P+1)$, all with the exact quadrature:

$$
M = B^\top W B, \qquad S = B^\top W B', \qquad K = B'^\top W B', \qquad T = \tilde B^\top W B,
$$

where $B_{qp} = \psi_p(\xi_q)$, $B'_{qp} = \psi_p'(\xi_q)$, $\tilde B_{qk} = \tilde\psi_k(\xi_q)$ (orthonormal Legendre) and $W = \operatorname{diag}(w_q)$.

## 2.2 Physical scaling

For element $e$ with length $h_e$, $x = x_e + (1+\xi)\, h_e/2$, so $dx = J\, d\xi$, $J = h_e/2$,
and $\dfrac{d}{dx} = \dfrac{1}{J}\dfrac{d}{d\xi}$. The operators used in the RHS become

| term | reference form | element form |
|---|---|---|
| mass $(\psi_i, u)$ | $M \hat u$ | $J M \hat u$ |
| convection $(\psi_i, \partial f/\partial x)$ | $S \hat f$ | $S \hat f$ (the $J$ and $1/J$ cancel) |
| viscous $(\psi_i', G)$ (after IBP) | $B'^\top W G_q$ | $B'^\top W G_q$ (same cancellation) |
| Laplacian $(\psi_i', u_x)$ | $K \hat u$ | $\frac{1}{J} K \hat u$ |
| SVV $(\psi_i', F * u_x)$ | $L \hat u$ | $\frac{1}{J} L \hat u$ |

## 2.3 Non-linear terms: projection and aliasing

$f(u) = u^2/2$ is evaluated at the $Q$ quadrature points and **projected** back onto
$\mathbb{P}_P$: $\hat f = M^{-1} B^\top W f_q$ (`Basis1D::project`). Then $(\psi_i, f_x) = S \hat f$. With
$Q = P+1$ the product $u^2$ (degree $2P$) is aliased; with $Q \ge (3P+3)/2$ (the "3/2
rule") the projection is exact.

Three algebraically equivalent, discretely different forms are implemented for
Burgers (`ConvectiveForm`): conservative $\left(v, (u^2/2)_x\right)$, non-conservative
$(v, u\, u_x)$, and skew-symmetric $\tfrac13 \left(v, (u^2)_x\right) + \tfrac13 (v, u\, u_x)$. The last one is the
"energy-stable" form for inexact quadrature. Their stability without SVV differs
enormously (see `docs/03`): this is a good reminder that SVV is a fix for the tail of
the spectrum, not for an aliasing-unstable discretisation.

For compressible NS, gradients of the *conserved* variables are computed exactly from
the modal expansion (they are polynomials), and $u_x$, $p_x$, $T_x$ follow by the
chain rule at the quadrature points. This avoids differentiating non-polynomial
quantities such as $u = m/\rho$.

## 2.4 Assembly and continuity

Global dof numbering: $\mathrm{dof}(e,p) = eP + p$. Element $e$'s right vertex $eP + P$
*is* element $(e+1)$'s left vertex $(e+1)P + 0$. That single line is the whole $C^0$
machinery in 1D. Periodic meshes take the index modulo $N_{el} P$.

Global matrices are assembled by scatter-adding element matrices; the mass matrix is
LU-factored once. Sizes here are at most a few hundred, so dense is fine (Nektar++ uses
static condensation and banded/iterative solvers for the same operation).

## 2.5 Boundary conditions

* **Dirichlet**: the rate of the boundary unknown is zero. Replace that row of $M$
  by a unit row and zero that entry of the residual. Exact and trivial.
* **Periodic**: only the dof map changes.
* Interface flux terms $[v\, G]$ from integration by parts are *not* added at element
  interfaces (standard $C^0$ Galerkin: the global test function is continuous, the weak
  form has only domain-boundary terms). They are dropped at the domain boundary too:
  irrelevant for Dirichlet rows, and consistent for periodic cases.

## 2.6 Time integration and the time step

Method of lines: the integrator sees $R(\hat u) = M_{bc}^{-1}\, r(\hat u)$. RK4 with

$$
\Delta t_{adv} = \mathrm{CFL}\, \frac{h_{min}}{|\lambda|_{max}}, \qquad
\Delta t_{diff} = 0.8\, \frac{2.785}{\rho(M^{-1} A)},
$$

where $h_{min}$ is the smallest physical GLL spacing and $A = \epsilon L_{svv}$, or $\nu K$ for physical viscosity. `burgers_svv` uses $|\lambda|_{max} = 1.2 \max\lvert u_0\rvert$ (20 % margin for the steepening profile); `cns1d_svv` uses $\max(|u| + c)$.

$\rho(M^{-1}A)$ is computed by power iteration (`CGAssembly::spectralRadius`), a cheap and
honest alternative to guessing $C\, P^4/h^2$. Both estimates are printed at start-up.

With the *split* SVV application the SVV term is implicit and drops out of the
restriction; the apps then only use the advective (and physical-viscous) limits.
The price is a first-order splitting error: the split result depends on $\Delta t$
(`docs/03` §3.3).

## 2.7 Initial conditions

Global $L^2$ projection: solve $M_{global}\, \hat u = \sum_e J_e B^\top W f(x_q)$. For a discontinuous
initial condition (Sod) this produces Gibbs oscillations at $t = 0$ which can already
make the pressure negative at $P = 8$; `--ic-smooth δ` replaces the jump by a
$\tanh\!\left((x - x_0)/\delta\right)$ profile. The projection error for smooth data is spectrally small:
the unit test projects a cubic exactly, and the entropy-wave case reaches $10^{-12}$.
