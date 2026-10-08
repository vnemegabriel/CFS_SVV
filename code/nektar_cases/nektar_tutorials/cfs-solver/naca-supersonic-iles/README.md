# Supersonic NACA implicit-LES — Nektar++ 5.9.0

M∞ = 1.5, spanwise-periodic slab, `NavierStokesCFE` with physical artificial
viscosity. Files:

| file | what it is |
|---|---|
| `make_naca_mesh.py` | Gmsh C-grid generator, extruded to a periodic slab |
| `session_naca_sup.xml` | Nektar++ session — every choice is commented in place |
| `run.sh` | mesh → NekMesh → solver → FieldConvert |

---

## 1. Two blockers in your current build

Both come from `/opt/nektar++/build/CMakeCache.txt`.

**`NEKTAR_USE_MPI:BOOL=OFF`.** This is the important one. Your solver has no
MPI. When you run `mpirun -np 4 CompressibleFlowSolver …`, MPI starts four
processes, but each one links Nektar's *serial* communicator, so each solves
the entire problem independently and writes over the others' output. It is not
a domain decomposition. The evidence was already in your earlier terminal
output: the session header printed four times. Under real MPI, only rank 0
prints it. A 3D LES is impossible until this is fixed.

**`NEKTAR_USE_MESHGEN:BOOL=OFF`.** NekMesh has no CAD kernel, so `.mcf`-driven
mesh generation is unavailable. Hence the Gmsh route below, which is the better
choice here anyway — you need explicit control of the wall-normal spacing.

Two lesser ones: `NEKTAR_USE_HDF5=OFF` means every rank writes its own `.fld`
file (unworkable above a few hundred ranks), and PETSc was built by ThirdParty
as MPI-uni, which is where the "you need full MPI version of PETSc" message
comes from.

### Rebuild

```bash
sudo apt install libopenmpi-dev openmpi-bin
cd /opt/nektar++/build
rm -rf ThirdParty/petsc*            # the MPI-uni PETSc must go
cmake -DNEKTAR_USE_MPI=ON \
      -DNEKTAR_USE_HDF5=ON  -DTHIRDPARTY_BUILD_HDF5=ON \
      -DNEKTAR_USE_SCOTCH=ON \
      -DNEKTAR_USE_PETSC=OFF \
      -DCMAKE_BUILD_TYPE=Release \
      ..
make -j"$(nproc)" install
```

PETSc is off because nothing here needs it — the explicit solver has no linear
system, and Nektar's own iterative solvers cover the implicit path. Leave it
off unless you specifically want PETSc preconditioners, in which case it must
be rebuilt against the same MPI. `Release` rather than `RelWithDebInfo` is
worth roughly 10–20% on a production run.

---

## 2. What "turbulent Navier–Stokes" can mean here

Nektar++ 5.9's `CompressibleFlowSolver` ships **no RANS model** — no
Spalart–Allmaras, no k–ω SST, no algebraic closure, and no wall model. I
checked the source tree; there is nothing to enable.

Turbulence in this solver is therefore **implicit LES**: you resolve the
energetic scales on a fine 3D mesh and let the DG scheme's numerical
dissipation act as the subgrid model. That is what this case is set up to do.
It has two consequences worth being explicit about in a thesis:

- **It must be 3D.** Two-dimensional "turbulence" has an inverse energy
  cascade and is a different physical system. The spanwise-periodic slab is the
  minimum honest geometry.
- **The mesh *is* the model.** With no explicit SGS term, resolution is not a
  convergence question, it is a modelling question. Section 4 quantifies it.

---

## 3. Physics and numerics, in brief

The reasoning for each setting lives in the comments of the two files. The
three decisions that matter most:

**Domain sizing is set by the Mach cone.** At M∞ = 1.5 the Mach angle is
μ = arcsin(1/M) = 41.8°. The leading-edge shock reaches the far-field height H
at x = H/tan μ downstream. You want it to leave through the *outlet*, where
`PressureOutflow` extrapolates cleanly at supersonic normal Mach, not through
the far field, where a Dirichlet free-stream state reflects it back onto the
aerofoil. Hence

> **H > L_x · tan μ**

The mesh script enforces this and refuses to build a domain that violates it.
Because it is satisfied, the plain Dirichlet far field is not a compromise —
at supersonic inflow all five characteristics enter, so the state is fully
determined and nothing is over-specified.

**Shock capturing has to be told what is a shock and what is turbulence.**
`ShockSensorType=Dilatation` keys on ∇·u, and `DucrosSensor=On` multiplies it
by θ²/(θ²+ω²), which → 1 in compression-dominated regions and → 0 in
vorticity-dominated ones. Without the Ducros filter, artificial viscosity is
sprayed across the turbulent boundary layer and wake and silently damps the
fluctuations the LES exists to resolve. This is the single most common way a
shock-capturing LES produces a plausible-looking but wrong answer.

`ShockCaptureType=Physical` requires a Navier–Stokes solver — the artificial
viscosity is added through the real diffusion operator. This is why the case
cannot be run as `EulerCFE`.

**Adiabatic wall, not isothermal.** At M = 1.5 the recovery temperature is
T_w/T∞ = 1 + r(γ−1)M²/2 = 1.40 with r = Pr^(1/3). An isothermal wall at T∞
would impose a 40% temperature deficit that does not physically exist.

---

## 4. Resolution and cost — read this before choosing Re

You asked for Re ≥ 10⁶. Here is what that costs, against Re = 6×10⁴ for
comparison. Wall units from the turbulent flat-plate correlation
c_f = 0.058 Re^−0.2, u_τ/U∞ = √(c_f/2), ℓ_ν = ν/u_τ with ν = 1/Re.

| | Re = 6×10⁴ | Re = 10⁶ |
|---|---|---|
| c_f | 6.4×10⁻³ | 3.7×10⁻³ |
| u_τ/U∞ | 0.057 | 0.043 |
| 1 wall unit, ℓ_ν/c | 2.9×10⁻⁴ | 2.3×10⁻⁵ |
| δ at trailing edge, /c | 0.041 | 0.023 |
| first element for y₁⁺ ≈ 1 (P=3) | 1.7×10⁻³ c | 1.4×10⁻⁴ c |
| **span 0.1c, in wall units** | **340** | **4270** |

That span row is the non-obvious one. A slab must be wide enough to decorrelate
across the span, which takes several near-wall streak spacings (≈100 wall units
each) — a working minimum is ~1000 wall units. At Re = 10⁶, 0.1c gives 4270 and
is comfortable. **At Re = 6×10⁴, 0.1c gives only 340 and is too narrow** — use
`--span 0.3`. Spanwise extent is not a free parameter; it scales with Re.

Now the cost, at P = 3 with element spacings Δx⁺ ≈ 160, Δz⁺ ≈ 80:

| | Re = 6×10⁴ | Re = 10⁶ |
|---|---|---|
| hexahedra | ~2×10⁴ | ~1.3×10⁶ |
| DOF (5 vars) | ~7×10⁶ | ~4×10⁸ |
| explicit Δt (c/U∞) | 7×10⁻⁵ | 6×10⁻⁶ |
| steps for 40 convective times | 5.5×10⁵ | 6.9×10⁶ |
| **estimated core-hours** | **~5×10²** | **~5×10⁵** |

Assumes ~2×10⁶ DOF-updates/core/s, which is optimistic-realistic for DG on
modern hardware; treat the last row as order-of-magnitude. The scaling is the
familiar wall-resolved-LES one, N ~ Re^1.86.

**So:** Re = 6×10⁴ is roughly a day on 32 cores. Re = 10⁶ is a national-facility
allocation — three orders of magnitude more, and that is before you discover
you need to repeat it at a second angle of attack. The `les` preset in the mesh
script is genuinely sized for Re ≈ 10⁶ (1.25×10⁶ hexes, y₁⁺ ≈ 0.45), so the
files are ready if you have the allocation. If you do not, `Re` in the session
file is a single line, and Re = 6×10⁴ with `--span 0.3 --level medium` gives a
defensible, publishable, reproducible result.

There is no third option that gets you Re = 10⁶ cheaply. Wall-modelled LES
would (N ~ Re rather than Re^1.86), but Nektar++ 5.9 has no wall model.

---

## 5. Transition — the case will not become turbulent on its own

A uniform initial condition contains no disturbance, and a supersonic boundary
layer is strongly stabilised by compressibility. Left alone this runs as a
steady laminar solution with a shock, and calling it an LES would be wrong.

Cheapest trip: broadband noise in the spanwise momentum, localised just aft of
the leading edge. `awgn(σ)` is Nektar's unit-variance white-noise function,
evaluated per quadrature point.

```xml
<E VAR="rhow" VALUE="rhoInf*(wInf + 0.02*awgn(1.0)*exp(-((x-0.1)^2 + y^2)/0.0025))" />
```

That seeds transition but convects away; it works because the resulting
turbulence is self-sustaining once established. If it is not — check after ~10
convective times — you need a permanent trip: a geometric step in the mesh near
x/c ≈ 0.05, or a steady wall-normal blowing/suction strip via `<FORCING>`.
State whichever you used; trip method is a reportable part of an LES setup.

---

## 6. Running it

```bash
chmod +x run.sh
./run.sh coarse 4      # ~34k elements, validates the whole pipeline
./run.sh medium 64
./run.sh les 512
```

Start with `coarse`. It exists to prove the mesh converts, the composites map
correctly, the periodic pair matches, and the boundary conditions do not blow
up — not to produce physics.

### Checklist before trusting any result

1. **Look at the artificial viscosity field.** `FieldConvert` writes it with the
   flow. It should be non-zero *only* in a thin band at the shock. If it is
   active in the boundary layer or wake, lower `mu0` and rerun — that is the
   LES being damped.
2. **Check the shock leaves through the outlet,** not the far field. Plot
   density gradient magnitude; if there is a reflected wave coming back off the
   top boundary, increase `--radius`.
3. **Check y₁⁺ from the actual solution,** not from the correlation above.
4. **Check the span decorrelates:** two-point spanwise correlation of u' should
   fall to ~0 by half the span. If it does not, the span is too narrow and the
   turbulence is artificially constrained.
5. **Average.** An instantaneous LES field means nothing on its own; the
   `AverageFields` filter is already configured.

### Note on your existing implicit tutorial case

Separately from all the above, `cfs-Euler-01/session_naca_im.xml` still has two
faults: it points at a restart file `naca_im_start.fld` that was never
generated, and its `<TIMEINTEGRATIONSCHEME>` is explicit SSP-RK3 where the
implicit solver needs `DIRKOrder2`. Both are described in the earlier
conversation; neither affects the files here.
