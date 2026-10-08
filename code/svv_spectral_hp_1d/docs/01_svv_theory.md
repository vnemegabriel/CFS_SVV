# 1. Spectral Vanishing Viscosity: from Tadmor to the C0 spectral/hp operator

This note follows Karniadakis & Sherwin (K&S), §6.5.2, pp. 348–352, and fills in the steps the book skips. Equation numbers refer to the book. Code references point to `src/03_SVV.h`.

## 1.1 The idea (Fourier, eq. 6.5.10–6.5.12)

Take the inviscid Burgers equation $u_t + \left(\tfrac{u^2}{2}\right)_x = 0$. A pure Galerkin/spectral discretisation has no dissipation, so once a shock forms the energy that should be dissipated at the discontinuity piles up in the highest modes and produces wiggles (and, often, blow-up). Adding an *artificial viscosity* $\epsilon\, u_{xx}$ cures the wiggles but destroys spectral accuracy: the error becomes $\mathcal{O}(\epsilon)$ everywhere, including where the solution is smooth.

Tadmor's (1989) observation: apply the viscosity **only to the high modes**,

$$
u_t + f(u)_x = \epsilon\, \frac{\partial}{\partial x}\left(F * u_x\right) \qquad (6.5.11)
$$

where $F*$ is a convolution that, in Fourier space, multiplies mode $k$ of $u_x$ by $\hat F_k$, with

$$
\hat F_k =
\begin{cases}
0, & |k| \le M_{SVV}, \\
1, & |k| > M_{SVV}.
\end{cases}
\qquad (6.5.12)
$$

* Modes below $M_{SVV}$ see **no viscosity at all**. For a smooth solution those are the only modes with energy, so the method keeps spectral accuracy.
* Modes above $M_{SVV}$ see a viscosity $\mathcal{O}(\epsilon)$. That is enough to remove the "spectral blocking" of energy at the grid scale and to make the scheme converge (in $L^2$) to the correct entropy solution.
* $\epsilon$ is a *vanishing* viscosity: theory uses $\epsilon \sim 1/P$, $M_{SVV} \sim 5\sqrt{P}$ (Maday, Ould Kaber & Tadmor 1993 for Legendre). The book's Burgers test uses $\epsilon = 1/16$, $M_{SVV} = 8$ with $P = 15$.

Exponential kernel used here (`svvKernel`, eq. 6.5.13):

$$
\hat F_k =
\begin{cases}
0, & k \le M_{SVV}, \\
\exp\!\left(-\dfrac{(k-P)^2}{(k-M_{SVV})^2}\right), & k > M_{SVV},
\end{cases}
\qquad \hat F_P = 1,\ C^\infty \text{ ramp.}
$$

The original 0–1 kernel (eq. 6.5.12) is available as `--kernel step` in the code; the sweep in `docs/03` shows it is worse (a sharp cut-off in mode space is itself a source of oscillation).

## 1.2 What "mode k" means in an element (the T matrix, p. 351)

In a spectral/hp element the natural basis is the **modified $C^0$ basis** (K&S 2.40):

$$
\psi_0 = \frac{1-\xi}{2}, \qquad
\psi_P = \frac{1+\xi}{2}, \qquad
\psi_p = \frac{1-\xi}{2}\,\frac{1+\xi}{2}\, P^{1,1}_{p-1}(\xi), \quad 0 < p < P.
$$

It is *semi-orthogonal*: the vertex modes overlap with everything. "Mode $p$" of this basis is not a frequency. The kernel needs a basis where mode $k$ = polynomial degree $k$: the orthonormal Legendre basis $\tilde\psi_k = \sqrt{\tfrac{2k+1}{2}}\, L_k$.
Both bases span $\mathbb{P}_P([-1,1])$, so there is an invertible change of basis `Basis1D::T`:

$$
\tilde u = T \hat u, \qquad T_{kp} = (\tilde\psi_k, \psi_p),
$$

because $\tilde u_k = (\tilde\psi_k, u) = \sum_p \hat u_p\, (\tilde\psi_k, \psi_p)$. Two facts checked by the unit tests:

* **Parseval**: $T^\top T = M$ (mass matrix), hence $T^{-1} = M^{-1} T^\top$. This is the identity the book states without proof on p. 351.
* The modal filter $\Phi = T^{-1} F T$, $F = \operatorname{diag}(\hat F_0, \dots, \hat F_P)$ (`modalFilter`), is *symmetric in the $L^2$ sense*: $(v, \Phi u) = (\Phi v, u)$, which is what makes the operators below symmetric.

>REVISAR
## 1.3 The discrete SVV operator (p. 351)

Weak form of the SVV term, after integration by parts in each element:

$$
\left(v,\ \epsilon\,(F * u_x)_x\right) \;\longrightarrow\; -\epsilon\,\left(v_x,\ F * u_x\right) \qquad (6.5.15)
$$

Per element, with the reference-element matrices

$$
M_{ij} = (\psi_i, \psi_j), \qquad S_{ij} = (\psi_i, \psi_j'), \qquad K_{ij} = (\psi_i', \psi_j'),
$$

the three steps are

1. Modal coefficients of $u_\xi$: $\hat w = M^{-1} S \hat u$ (since $(\psi_i, u_\xi) = \sum_j S_{ij} \hat u_j$).
2. Filter in the orthonormal space: $\hat w_f = T^{-1} F T \hat w$.
3. Inner product with $v_\xi$: $(\psi_i', w_f) = \sum_j S_{ji}\, \hat w_{f,j} = (S^\top \hat w_f)_i$.

Hence

$$
L_{svv} = S^\top T^{-1} F T M^{-1} S = S^\top M^{-1} T^\top F T M^{-1} S \qquad (\texttt{svvOperator}, \text{ form "book"})
$$

Properties (all verified in `tests/test_svv_operator.cpp`):

* symmetric, positive semi-definite;
* $F = I$ gives exactly the Laplacian $K$ (because $S^\top M^{-1} S = K$: the projection of $u_\xi$ onto $\mathbb{P}_P$ is exact);
* **null space** = polynomials of degree $\le M_{SVV} + 1$ (their derivative has degree $\le M_{SVV}$, killed by the kernel). SVV is literally invisible to resolved modes.
* $u^\top L_{svv}\, u \le u^\top K u$: it dissipates less than the Laplacian.

Physical scaling: on an element of length $h$ each $d/dx$ brings a $2/h$ and the integral brings $h/2$, so

$$
(v_x, F u_x)_{\text{element}} = \frac{2}{h}\, L_{svv,\text{ref}} = \frac{1}{J}\, L_{svv,\text{ref}}.
$$

Kirby's variant (eq. 6.5.17) filters the solution before and after a Laplacian: $L_2 = \Phi^\top K \Phi$ (form "kirby"). Its null space is degree $\le M_{SVV}$, it is also symmetric positive semi-definite, and it is less dissipative in the sweep of `docs/03`.

## 1.4 A structural fact with big consequences: vertex rows are zero

The derivative of a vertex mode is a constant = Legendre mode 0. Every kernel has $\hat F_0 = 0$. Therefore

$$
(L_{svv})_{0,j} = (L_{svv})_{P,j} = 0,
$$

and, by symmetry, the vertex columns too (unit test `vertex rows and columns of L_svv are zero`).

**SVV never produces a residual on a vertex degree of freedom, and never reads one.**
In a $C^0$ discretisation the vertex dofs are the only ones shared between elements, so the SVV operator, by itself, does not couple elements at all. Two consequences:

1. Continuity is never at risk: a purely element-local application of SVV is possible and keeps the $C^0$ solution $C^0$.
2. If the operator is used as a weak-form residual through the **consistent global mass matrix** ($M\, \tfrac{du}{dt} = r_{conv} - \epsilon L u$), the bubble-mode residuals of element $e$ still move the shared vertex, because $M^{-1}$ couples bubbles and vertices. The neighbour must then absorb a moving vertex value with a function that is $L^2$-orthogonal to its own bubble space, which is a spiky, high-degree "lift". In a problem without physical viscosity nothing damps that lift. This makes the interface spikes documented in `docs/03` worse in the Galerkin application than in the split one at equal $\epsilon$.

The code therefore offers two applications (`SVVApplication`):

* **Galerkin**: the formula of p. 351 as a residual (the literal reading of the text);
* **Split**: after each explicit step, in each element take the bubble rows of the implicit-Euler step of the *physical* element system, $(J M)(u^{n+1} - u^n) = -\Delta t\, \tfrac{\epsilon}{J}\, L\, u^{n+1}$, with the vertex frozen, and divide by $J$:

  $$
  \left(M_{bb} + \Delta t\, \frac{\epsilon}{J^2}\, L_{bb}\right) u_b^{n+1} = M_{bb}\, u_b^{n}
  $$

  Because $L$ has no vertex coupling this is the exact bubble block of an implicit-Euler SVV step. It is unconditionally stable and has no interface coupling. It is a first-order splitting, so its result depends on $\Delta t$ (`docs/03` §3.3).

  > **Correction (September 2026).** The code used $\epsilon/J$ in place of $\epsilon/J^2$, which applies $\epsilon J$ instead of $\epsilon$. The operator tests could not see it, because they never compared the filter with the physical element system; the test `split SVV step = implicit Euler of the physical element system` now does.

The interface coupling through $M^{-1}$ explains why the Galerkin application is worse than the split one at equal $\epsilon$. It does not explain everything: at the book's amplitude the split application also shows interface spikes and high-mode growth in the elements next to the shock (`docs/03` §3.4).

In an incompressible solver (Nektar++ `VelocityCorrectionScheme`) SVV is added to a global implicit Helmholtz operator *together with the physical Laplacian*, which does act on the vertices and damps the lift. That is why the issue is invisible there and why it matters for an explicit, inviscid or nearly inviscid compressible solver.

## 1.5 Time-step consequences

The SVV term is a filtered Laplacian, so for an explicit scheme it adds a diffusive restriction. `CGAssembly::spectralRadius` computes $\rho(M^{-1} L)$ by power iteration and the apps use

$$
\Delta t \le 0.8 \cdot \frac{2.785}{\epsilon\, \rho}
$$

(RK4's negative real-axis limit). For the Burgers test that is $\Delta t \approx 10^{-3}$, versus $2.5 \times 10^{-3}$ from the advective CFL: SVV in the Galerkin form more than halves the time step. The split form is implicit and free of this restriction, which is another practical argument for it.

## 1.6 Reading list

* Tadmor, *Convergence of spectral methods for nonlinear conservation laws*, SINUM 1989.
* Maday, Ould Kaber, Tadmor, *Legendre pseudospectral viscosity method for nonlinear conservation laws*, SINUM 1993 (kernel 6.5.13, $\epsilon \sim 1/P$, $M \sim 5\sqrt{P}$).
* Karamanos & Karniadakis, *A spectral vanishing viscosity method for LES*, JCP 2000.
* Kirby & Karniadakis, *Coarse resolution turbulence simulations with SVV-LES*, J. Fluids Eng. 2002 (the source of Fig. 6.27).
* Nektar++: `StdRegions/StdExpansion::v_SVVLaplacianFilter`, and `SetUpSVV()` in `IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.cpp`.
