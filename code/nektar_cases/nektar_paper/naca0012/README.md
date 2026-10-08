# NACA0012 — Nektar++ CompressibleFlowSolver

2D compressible flow over NACA0012. Two configurations share one mesh:

| Config | File | Physics |
|---|---|---|
| Viscous | `conditions.xml` | `NavierStokesCFE`, laminar, Re_c = 500, α = 1.25° |
| Inviscid | `conditions-euler.xml` | `EulerCFE`, α = 1.5° — reproduces Fig. 4(b) of the Nektar++ 2015 CPC paper |

Nektar++ 5.9.0.

---

## Files

| File | Contents |
|---|---|
| `naca0012.xml` | Geometry only: vertices, edges, elements, curved edges, composites. No `<EXPANSIONS>`. |
| `conditions.xml` | Viscous case: expansions, parameters, solver info, BCs, filters. |
| `conditions-euler.xml` | Inviscid case, same structure. |
| `nekcase.py` | Parses/evaluates session parameters, derives freestream reference. Imported by the others. |
| `dt_estimate.py` | Prints a `TimeStep` and `NumSteps` for a given mesh + state. |
| `viz.sh` | Checkpoints → VTU with primitive variables, tagged `ns_`/`euler_`. |
| `forces.py` | `forces.fce` → Cl/Cd + `forces.png`. |

Both XMLs are required on every command line — the expansion definition lives in the conditions file, not the mesh.

Composites: `C[1]` = 3002 quads (domain), `C[5900]` = 77 aerofoil edges, `C[5901]` = 113 farfield edges.

---

## Running

```bash
mpirun -n 4 CompressibleFlowSolver naca0012.xml conditions.xml         # viscous
mpirun -n 4 CompressibleFlowSolver naca0012.xml conditions-euler.xml   # inviscid
```

Changing the rank count: delete `naca0012_xml/` first, or stale partitions from the previous count are left behind.

Any file named `naca0012.opt` is auto-loaded by session name without appearing in the banner. `_0.chk/Info.xml` lists every file the session actually read.

---

## Parameters

### State — the only three knobs

| Parameter | Effect |
|---|---|
| `Mach` | Sets `uInf`, `vInf`, `Uref`. `pInf` and `rhoInf` are fixed, so `cInf` = 340.29 m/s always and Mach scales velocity directly. Raising it tightens the convective timestep limit. |
| `alpha` | Angle of attack, degrees. Splits `Uref` into `uInf`/`vInf`. Also used by `forces.py` to rotate forces into wind axes. |
| `Reynolds` | Sets `mu = rhoInf*Uref*Lref/Reynolds`. Raising it *relaxes* the viscous timestep limit. Ignored by `EulerCFE`. |

Everything else is derived:

```xml
<P> TInf  = pInf/(GasConstant*rhoInf)     </P>   <!-- 288.14 K   -->
<P> cInf  = sqrt(Gamma*GasConstant*TInf)  </P>   <!-- 340.29 m/s -->
<P> uInf  = Mach*cInf*cos(alpha*PI/180)   </P>
<P> vInf  = Mach*cInf*sin(alpha*PI/180)   </P>
<P> mu    = rhoInf*Uref*Lref/Reynolds     </P>
```

The case is fully dimensional (real sea-level air). Reynolds number is imposed by choosing an artificial `mu`. Do not set `thermalConductivity` and `Pr` together — that is a fatal error.

### Fixed physical constants

| Parameter | Value | Effect |
|---|---|---|
| `Gamma` | 1.4 | Ratio of specific heats. |
| `GasConstant` | 287.058 | Sets `Cp = Gamma/(Gamma-1)*R`, `Cv = Cp/Gamma`. |
| `pInf` | 101325 | Reference pressure; also the Cp offset. |
| `rhoInf` | 1.225 | Reference density. |
| `Pr` | 0.72 | Thermal conductivity = `Cp*mu/Pr`. |
| `Lref` | 1.0 | Chord (measured 0.9992). |
| `Twall` | — | Wall temperature. Only read by `WallViscous`. Default 300.15. |

### SOLVERINFO

| Property | Value used | Other valid values | Effect |
|---|---|---|---|
| `EQType` | `NavierStokesCFE` / `EulerCFE` | `NavierStokesImplicitCFE` | Viscous terms on/off. Implicit variant removes the viscous timestep limit. |
| `DiffusionType` | `InteriorPenalty` | `LDGNS` | Viscous discretisation. Must be one of these two — plain `LDG` is a different class and segfaults. Omit for `EulerCFE`. |
| `AdvectionType` | `WeakDG` | — | |
| `UpwindType` | `HLLC` | `Roe`, `ExactToro`, `AUSM0-3`, `HLL`, `LaxFriedrichs` | Riemann solver. |
| `ShockCaptureType` | `Physical` / `NonSmooth` | `Off` | `NavierStokesCFE` accepts all three. `EulerCFE` accepts only `NonSmooth` and `Off`. |
| `ShockSensorType` | `Dilatation` | `Modal` | Only read when `ShockCaptureType=Physical`. |
| `DucrosSensor` | `On` | `Off` | Suppresses AV in vortical regions. `Physical` only. |
| `ViscosityType` | `Constant` | `Variable` | `Variable` enables Sutherland's law. |
| `Projection` | `DisContinuous` | — | |
| `LocalTimeStep` | not used | `True` | Per-element timestep for steady state. Requires `TimeStep = 0` and a non-zero `CFL`. Diverged here at CFL 0.2 → 0.02. |

### Shock capturing parameters

| Parameter | Read by | Effect |
|---|---|---|
| `mu0` | both | Scaling constant. `mu_av ≈ mu0*(h/P)*rho*(U+c)*sensor`, peaking near 5.6 Pa·s on this mesh at `mu0=1.0`. Raise for more damping, lower for a larger timestep. |
| `SensorOffset` | `NonSmooth` | Sensor threshold offset. Default 1. |
| `Skappa`, `Kappa` | `Physical` + `ShockSensorType=Modal` | Unused with the default `Dilatation` sensor. |

Artificial viscosity adds to `mu` and therefore tightens the viscous timestep limit. At Re = 500 (`mu` = 0.667) that is a ~4× bump and survivable; at Re = 5000 (`mu` = 0.0667) it is ~80× and the run diverges.

### Boundary conditions

Region 0 = aerofoil, region 1 = farfield (15 chords).

| `USERDEFINEDTYPE` | Effect |
|---|---|
| `WallAdiabatic` | No-slip, zero heat flux. Used for the viscous case. |
| `WallViscous` | No-slip, isothermal at `Twall`. |
| `Wall` | Slip / symmetry. Used for the Euler case. On a viscous run it produces no boundary layer. |
| `RiemannInvariant` | Characteristic farfield, non-reflecting. Used for both. |
| `PressureOutflow` | Static pressure imposed. |
| none | Hard Dirichlet. Reflects acoustic waves. |

### Filters

`Checkpoint` writes `.chk` at `OutputFrequency`. `AeroForces` on `B[0]` writes `forces.fce`. Use these rather than `IO_CheckSteps`, which is deprecated.

---

## Timestep

Explicit RK4 must satisfy a convective and a diffusion limit simultaneously:

```
dt_conv = 1.0   * h / ( (U+c) * (2P+1) )
dt_visc = 0.091 * h² / ( nu * (2P+1)² )

1/dt_max = 1/dt_conv + 1/dt_visc
TimeStep = 0.5 * dt_max
```

`h` = 1.171e-2 (minimum edge on this mesh), `P` = NUMMODES − 1 = 3. `Cv` = 0.091 was calibrated on this case; recalibrate if the mesh changes.

```bash
python3 dt_estimate.py                              # from conditions.xml
python3 dt_estimate.py --mach 1.5 --re 1000 --p 4   # override
python3 nekcase.py                                  # show derived freestream
```

At Re = 500 the viscous limit (4.67e-7) is 6× more restrictive than the convective limit (2.73e-6). Sizing from CFL alone gives NaN within ~100 steps.

Measured on this mesh, P = 3, M = 0.8:

| Config | Re | dt | Shock capturing | Result |
|---|---|---|---|---|
| NS | 500 | 3e-7 | Off | clean |
| NS | 500 | 5e-7 | Off | NaN |
| NS | 500 | 2e-7 | Physical | clean (10 000 steps) |
| NS | 5000 | 2e-6 | Off | clean |
| NS | 5000 | 1e-6 | Physical | NaN |
| Euler | — | 4e-7 | NonSmooth | NaN |
| Euler | — | 2e-7 | NonSmooth | clean |

Settled values: **2e-7** for both configurations.

One convective time is `c/U` = 3.673e-3 s at M = 0.8. A converged run needs 20–30 of them: ~370 000–460 000 steps. Measured throughput is 60 ms/step on 4 ranks, so 6–8 hours.

---

## Mesh envelope

Measured: first wall-normal cell 2.24 % of chord (1.92e-2 to 2.78e-2), surface spacing 1.17e-2 to 3.13e-2 over 77 elements, farfield at 15 chords. Near-isotropic — no boundary-layer clustering.

Against δ₉₉/c ≈ 5/√Re:

| Re | δ₉₉/c | Elements across BL |
|---|---|---|
| 500 | 0.224 | 10.0 |
| 1 000 | 0.158 | 7.1 |
| 5 000 | 0.071 | 3.2 (needs P ≥ 4) |
| 1e5 | 0.016 | 0.7 |
| 1.9e7 (real air) | 0.001 | 0.05 |

**Re ≈ 5000 is the ceiling.** Higher requires a NekMesh regeneration with a stretched wall layer.

---

## Post-processing

### Checkpoints → VTU

```bash
bash viz.sh                                # viscous case, every checkpoint
bash viz.sh 5                              # only checkpoint 5
COND=conditions-euler.xml bash viz.sh      # Euler case
KEEP=1 bash viz.sh                         # do not delete previous output
```

Output goes to `vtk/`, tagged by solver — `ns_*` or `euler_*` — read from `EQType` in the
conditions file. Each run deletes **its own** previous files first, so stale checkpoints
never linger, while the other solver's results are left untouched. Numbers are zero-padded
so the ParaView series stays ordered.

| Output | Contents |
|---|---|
| `vtk/<tag>_0005.vtu` | Volume field: `u v p T Cp Mach` + conservative |
| `vtk/<tag>_surf_0005_b0.vtu` | Aerofoil surface, conservative only |
| `vtk/<tag>_INFO.txt` | Solver, conditions file, Mach, alpha, Re, `pInf`, `q`, timestamp |

Under the hood:

```bash
FieldConvert -e naca0012.xml conditions.xml naca0012_5.chk out.vtu
```

`-e` writes equispaced points. Without it each P = 3 element is drawn as a single linear cell.

Checkpoints store only `rho, rhou, rhov, E`. `viz.sh` derives the rest with chained `fieldfromstring` modules:

| Field | Expression |
|---|---|
| `u`, `v` | `rhou/rho`, `rhov/rho` |
| `p` | `(Gamma-1)*(E - 0.5*(rhou²+rhov²)/rho)` |
| `Cp` | `(p - pInf)/q`, with `q = 0.5*rho*U²*c` |
| `T` | `p/(R*rho)` |
| `Mach` | `sqrt(u²+v²)/sqrt(Gamma*p/rho)` |

Constants come from `conditions.xml` via `nekcase.py`, so they follow the Mach number automatically.

### ParaView

Open `vtk/ns_....vtu` or `vtk/euler_....vtu` — the dots load the set as an animatable time series. Colour by `Mach` or `Cp`.

For the shock: **Contour** on `Mach` at 1.0.

### Surface Cp

`viz.sh` also writes `<tag>_surf_NNNN_b0.vtu` (boundary 0). The `extract` module discards derived fields, so these carry conservative variables only. Rebuild Cp with a ParaView **Calculator**:

```
(0.4*(E-0.5*(rhou^2+rhov^2)/rho) - pInf) / q
```

`viz.sh` prints the correct `pInf` and `q` for the current state at the end of every run. Then **Plot Data** against `Points_X`.

### Forces

```bash
python3 forces.py                      # forces.fce + conditions.xml
python3 forces.py run2.fce cond2.xml
```

`AeroForces` writes body-axis components. Lift and drag are those rotated by α:

```
Cd = ( Fx*cos(alpha) + Fy*sin(alpha)) / q
Cl = (-Fx*sin(alpha) + Fy*cos(alpha)) / q       q = 0.5*rho*U²*c
```

Prints mean and standard deviation over the second half of the history, and warns below ten convective times. Writes `forces.png`.

---

## Reproducing Fig. 4(b)

The paper figure is **Euler**, Ma∞ = 0.8, α = 1.5°, peak Mach 1.37, showing a strong shock on the upper surface and a weak one on the lower.

Use `conditions-euler.xml`. It cannot be reproduced with `conditions.xml`: at Re = 500, ν = 0.545 m²/s and the boundary layer is 22 % of chord, so the flow never reaches M = 1 anywhere. Measured peak Mach in the viscous case is 0.92.

Progress of the Euler run (peak Mach in the field):

| Step | Convective times | Mach max | M > 0.95 | M > 1 |
|---|---|---|---|---|
| 20 000 | 1.09 | 0.9613 | 758 | 0 |
| 25 000 | 1.36 | 0.9870 | 2 147 | 0 |

Monotonically approaching sonic and the near-sonic region is growing fast (758 → 2147 points in 0.27 convective times). The supersonic pocket has not opened yet at 1.36 convective times; budget 20–30 for a steady solution.

Measured cost: 60 ms/step on 4 ranks, so 370 000 steps ≈ 6 hours.
