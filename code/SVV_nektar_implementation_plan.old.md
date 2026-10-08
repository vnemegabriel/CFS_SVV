# SVV in Nektar++ CompressibleFlowSolver — Implementation Plan

*Thesis scope: implement a Spectral Vanishing Viscosity (SVV) stabilizer in Nektar++'s compressible Navier–Stokes solver and benchmark it against the library's existing stabilization mechanisms on Taylor–Green vortex, Shu–Osher, and forward-facing step across sub- and supersonic regimes.*

Survey basis: `TESIS/code/ndg_methods` (Hesthaven–Warburton Codes1.1, MATLAB), `TESIS/code/ESDG` (Julia, entropy-stable DG + IDP for CNS), `TESIS/code/3d_mortensen_dns` (pseudo-spectral TGV DNS reference), Nektar++ **5.9.0** source at `/opt/nektar++` (WSL).

---

## 1. What already exists (survey results)

### 1.1 Nektar++ compressible solver architecture

`solvers/CompressibleFlowSolver/` is pure DG (`eDiscontinuous` enforced), explicit RK or implicit. The RHS assembly in `CompressibleFlowSystem::DoOdeRhs` is:

```
outarray = −Advection(WeakDG + Riemann solver)
         + Diffusion(LDGNS or InteriorPenalty)   [NavierStokesCFE::v_DoDiffusion]
         + Forcing
```

Four stabilization mechanisms already present — these are your **comparison baselines**, all usable with zero code changes:

| Mechanism | Session key | Where it acts |
|---|---|---|
| Over-integration dealiasing | `SPECTRALHPDEALIASING` | `GetFluxVectorDeAlias` |
| Exponential modal filter | `ExponentialFiltering` + `FilterAlpha/Exponent/Cutoff` | applied to solution every step in `DoOdeProjection` via `ExpList::ExponentialFilter` |
| Laplacian artificial viscosity (Persson–Peraire sensor) | `ShockCaptureType = NonSmooth` | `ArtificialDiffusion` factory → own LDG operator on conserved vars |
| Physical AV (div/curl sensor augments µ) | `ShockCaptureType = Physical` (IP diffusion) | `m_varConv->SetAv` inside `v_DoDiffusion` |

### 1.2 SVV machinery already in the library

`StdRegions` implements per-element SVV kernels via orthogonal modal projection (`StdQuadExp::v_SVVLaplacianFilter` and siblings for Tri/Hex/Tet/Prism/Pyr, using `eOrtho_A/B` bases — the collapsed-coordinate problem on triangles is already solved). Three kernels selectable through `StdMatrixKey` const factors:

- `eFactorSVVCutoffRatio` + `eFactorSVVDiffCoeff` — classical exponential kernel (Maday–Tadmor / Karniadakis–Kirby)
- `eFactorSVVPowerKerDiffCoeff` — Moura's power kernel
- `eFactorSVVDGKerDiffCoeff` — Moura/Mansoor DG-mimicking kernel

**The gap (your thesis contribution):** these kernels are only reachable through CG Helmholtz/Laplacian matrix operators used by the *incompressible* `VelocityCorrectionScheme` implicit viscous solve. Nothing connects them to the explicit DG compressible path.

### 1.3 The insertion point

`ArtificialDiffusion` (base class, `CompressibleFlowSolver/ArtificialDiffusion/`) is a small, clean plugin: it owns an LDG diffusion operator over the conserved variables and supplies a flux vector

```cpp
viscousTensor[j][i] = mu(x) * qfield[j][i]      // GetFluxVector
```

where `qfield = ∇u` comes from the LDG auxiliary solve. `NonSmoothShockCapture` merely defines the scalar `mu(x)` from the Persson sensor. **SVV drops in by replacing the pointwise multiply with a per-element modal operation on `qfield`.**

### 1.4 Local repos — role in the plan

- **ndg_methods** — prototyping sandbox. `Codes1D/Filter1D.m` is the exponential cutoff filter (NDG book §5.3, exactly what you're studying); `CFD1D/Euler*` gives 1D Euler with slope limiting (Sod, Shu–Osher); `CFD2D/EulerShock2D.m` + `ForwardStepBC2D/IC2D.m` is the Mach 3 forward step; `CurvedCNS2D.m` is compressible NS. SVV prototypes go here first.
- **ESDG** — reference for entropy-consistent formulation and hard test cases: `examples/IDP/dg1D_euler_shuosher.jl`, 2D CNS shocktube/double-Mach with IDP limiting. Guides the (optional, phase-2) entropy-variable variant of SVV and provides independent reference solutions.
- **3d_mortensen_dns** — pseudo-spectral TGV DNS (`mortensens.py`, `tgv_final.h5`, spectrum plot): reference for TGV kinetic-energy dissipation rate and spectra.

---

## 2. Formulation

Add to the semi-discrete conserved system a spectral viscosity term

$$\partial_t u + \nabla\!\cdot F(u) = \nabla\!\cdot\big(\varepsilon\, Q_N * \nabla u\big),\qquad \varepsilon = \mu_{SVV}\,\frac{\bar c\, h}{N}$$

where per element the kernel acts modally: expand each component of $\nabla u$ in the orthogonal (Dubiner/tensor) basis, multiply mode $k$ by $\hat Q_k$ with $\hat Q_k = 0$ for $k \le P = r_{cut} N$ and $\hat Q_k \to 1$ smoothly for $k \to N$ (exponential kernel; power/DG kernels as alternatives). This preserves design-order accuracy on smooth solutions (kernel vanishes on resolved modes) while dissipating the marginally-resolved tail — the classical Tadmor result, in the element-local form used by Nektar++'s own incompressible SVV.

Practical DG realization (identical structure to `ArtificialDiffusion::GetFluxVector`):

1. LDG auxiliary solve gives `qfield[j][i]` (already in place).
2. Per element: forward-transform each `qfield` component to the orthogonal basis, scale coefficients by $\varepsilon\,\hat Q_k$, backward-transform. Mirror the code in `v_SVVLaplacianFilter` (which does exactly this projection) but as a standalone filter of the gradient rather than inside a matrix op.
3. Return the filtered, scaled tensor as the viscous flux; LDG divergence + traces complete the term.

Design choices to characterize (this is thesis material, not overhead):

- **Kernel**: exponential vs power vs DG kernel — all three already coded in `StdRegions`, so the comparison is nearly free.
- **Amplitude scaling**: $\varepsilon \propto h/N$ vs $h/N^2$; constant vs velocity-scaled ($\bar c = |u|+a$ locally).
- **Variables acted on**: conserved $u$ (baseline) vs entropy variables $v(u)$ (entropy-consistent SVV, guided by ESDG) vs primitive.
- **Shock handling**: SVV is an $O(h/N)$ viscosity on high modes only — it stabilizes under-resolution/turbulence but is *not* sufficient at strong shocks (FFS Mach 3). Plan for sensor-gated amplification: $\varepsilon \to \varepsilon\,(1 + \beta\, s_e)$ with the Persson sensor $s_e$ already available from `NonSmoothShockCapture`. Be explicit about this division of labor in the thesis.

---

## 3. Implementation phases

### Phase 0 — MATLAB prototype in ndg_methods (2–3 weeks, parallel to NDG ch. 5 study)

1. 1D: add `SVVFilter1D.m` (modal kernel on $\partial_x u$, structure of `Filter1D.m`) and an `EulerRHS1D` variant with the SVV term via the existing LDG-style gradient. Cases: smooth advection (convergence — order preserved?), Burgers (spectra), Sod, Shu–Osher.
2. Sweep kernel type, $r_{cut} \in \{0.5, 0.7, 0.9\}$, $\mu_{SVV}$; establish the $\varepsilon(N,h)$ scaling that survives Shu–Osher without a limiter (or identify where sensor gating becomes mandatory).
3. 2D check: bolt the same construction onto `EulerShock2D` forward step.

Deliverable: validated kernel + scaling law + failure-mode map. This de-risks every downstream C++ decision and produces thesis figures directly.

### Phase 1 — Nektar++ groundwork (1–2 weeks)

1. Build 5.9.0 in WSL (`-DNEKTAR_SOLVER_COMPRESSIBLE_FLOW=ON`; pin to the release tag, work on a branch — 5.9 has active ALE churn in `CompressibleFlowSystem`).
2. Run baselines: isentropic vortex (convergence harness), subsonic TGV, FFS with `NonSmooth`, one `ExponentialFiltering` case. Learn the session-XML surface you'll extend.

### Phase 2 — `SVVDiffusion` plugin (3–4 weeks)

Touch list:

| File | Change |
|---|---|
| `ArtificialDiffusion/SVVDiffusion.{h,cpp}` | **new** — registered in `ArtificialDiffusionFactory` as `"SVV"`; overrides `GetFluxVector` with per-element modal kernel on `qfield`; reads `SVVDiffCoeff`, `SVVCutoffRatio`, kernel selector from session |
| `CompressibleFlowSolver/CMakeLists.txt` | add sources |
| `EquationSystems/CompressibleFlowSystem.cpp` | instantiate SVV via its own SolverInfo key (`SpectralVanishingViscosity = True`) rather than through `ShockCaptureType`, so SVV can coexist with — or replace — shock capture; keep `ShockCaptureType = SVV` as an alias if convenient |
| `EulerCFE.h` / `NavierStokesCFE.h` | extend `v_SupportsShockCaptType` if using the alias route |

Implementation notes:

- Do the modal filtering element-wise inside the plugin (loop `m_fields[0]->GetExp(n)`, orthogonal `FwdTrans`/scale/`BwdTrans` of the gradient), copying the projection pattern from `StdQuadExp::v_SVVLaplacianFilter`. Quads + tris first; hex for 3D TGV.
- Provide the `DiffuseCoeffs` path too (mirrors `v_DoArtificialDiffusionCoeff`) so implicit/ALE paths don't break.
- Sensor gating: reuse the Persson sensor from `NonSmoothShockCapture` (either composition or a small shared helper).
- Time-step: the added term is second-order — explicit `dt` scales like $h^2/(N^4 \varepsilon)$. With $\varepsilon \sim h/N$ this stays subdominant to the advective CFL at moderate N; verify in `GetElmtTimeStep` and document.
- Unit/regression tests in the solver's test format (isentropic vortex with SVV on: error must match design order).

### Phase 3 — Validation & comparison campaign (4–6 weeks)

Matrix: SVV vs the four existing baselines (§1.1), across Mach.

| Case | Regime | Metrics | Reference |
|---|---|---|---|
| Isentropic vortex | smooth, Ma 0.5 | $L_2$ convergence, $p = 2..8$ | analytic |
| TGV Re 1600 | Ma 0.1 (+ 0.5, 1.25 if time) | $-dE_k/dt$, enstrophy, energy spectra; DOF study $p = 3..7$ | `mortensens.py` DNS + literature (DeBonis, van Rees) |
| Shu–Osher | Ma 3 shock–entropy | density profile vs converged reference; oscillation amplitude | ESDG `dg1D_euler_shuosher.jl` / WENO reference (run on pseudo-1D strip mesh if 1D CFS support is limited) |
| Forward-facing step | Ma 3 | robustness (survives to $t=4$?), density/vorticity fields, entropy production | ndg_methods `EulerShock2D` + Woodward–Colella |

Cross-cutting measurements: CPU time per DOF per step (SVV's extra cost = one ortho transform pair per gradient component), maximum stable CFL, minimum resolution that survives each case, sensitivity to $\mu_{SVV}$ and $r_{cut}$.

### Phase 4 — (optional, if schedule allows) Entropy-variable SVV

Apply the kernel to $\nabla v(u)$ with entropy variables from ESDG's `EntropyStableEuler` formulas; demonstrates guaranteed entropy dissipation and is the clearest novelty claim beyond "port of SVV to compressible DG". Cleanly separable — do not let it block Phase 3.

---

## 4. Risks and mitigations

1. **SVV insufficient at strong shocks** — expected; sensor gating is in-plan, and the honest framing (SVV = under-resolution stabilizer, sensor = shock handler) strengthens rather than weakens the thesis.
2. **dt collapse from the viscous term** — controlled by $\varepsilon \sim h/N$ scaling; measured explicitly in Phase 3.
3. **5.9.0 API churn (ALE, `m_meshDistorted` branches)** — pin the tag; keep the entire change as one reviewable diff (also serves as a thesis appendix and a potential upstream MR).
4. **1D support in CompressibleFlowSolver** — Shu–Osher falls back to a thin 2D strip with periodic lateral BCs; prototype fidelity is covered by ndg_methods/ESDG anyway.
5. **Scope creep toward the full hybrid library** — this plan deliberately produces a *Nektar++-native* plugin, not a Python library; the earlier library architecture becomes future work in the thesis narrative.

## 5. Indicative schedule (from Aug 2026)

| Month | Milestone |
|---|---|
| Aug | Finish NDG ch. 5; Phase 0 1D prototype + parameter sweeps |
| Sep | Phase 0 2D forward step; Phase 1 Nektar++ builds + baselines |
| Oct | Phase 2: `SVVDiffusion` implemented, vortex regression passing |
| Nov | Phase 3: TGV + Shu–Osher campaigns |
| Dec | Phase 3: FFS + Mach sweep; consolidated comparison figures |
| Jan | Buffer / Phase 4 / writing |
