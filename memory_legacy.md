**Purpose & context**

Valen is a final-year mechanical engineering student with a minor in computational mechanics, working on a doctoral/engineering thesis centered on implementing and characterizing a **Spectral Vanishing Viscosity (SVV)** stabilization operator within Nektar++'s `CompressibleFlowSolver` for transonic regimes. The thesis contribution is framed as both *implementation* (SVV is entirely absent from `CompressibleFlowSolver` in the current codebase) and *characterization*, benchmarked against a second-order finite volume solver (OpenFOAM) on a transonic NACA 0012 case, with a planned computational campaign on the ARCHER2 HPC cluster.

The broader motivation is democratizing high-accuracy CFD for resource-constrained researchers, with an interest in hybrid library design combining custom spectral code with existing frameworks (Nektar++, Nek5000).

**Core references established:**
- Karniadakis & Sherwin (2nd ed., 2005) — primary reference for discretization, aliasing, SVV lineage (SVV: §6.5.2; compressible DG: §10.5.1)
- Kirby & Sherwin (2006, CMAME), Moura et al. (2016), Persson & Peraire (2006), Mateo-Gabín et al. (2022), Manzanero et al. (2020) — SVV in compressible/DG contexts
- Tadmor (1989 SIAM + CIME notes) — entropy-theoretic SVV foundations
- Wang et al. (2013) — work-precision diagram methodology
- Vermeire et al. (2017) — methodological precedent
- Slotnick (2014) — CFD vision/context

**Current state**

- Nektar++ 5.10.0 installed from source on Ubuntu 22.04 (WSL2), at `~/nektar-install`, with `lib64` as the library path. Build includes `CompressibleFlowSolver`, `IncNavierStokesSolver` (SVV-enabled), and `ADRSolver` (SVV test bed for scalar problems).
- `ctest` baseline run completed; parallel test failures from OpenMPI/SMT interaction resolved via `OMPI_MCA_hwloc_base_use_hwthreads_as_cpus=1` persisted to `~/.bashrc`.
- SVV grep verification confirmed: SVV kernels (power, DG, classical Tadmor) exist in `library/StdRegions` and are exposed through `IncNavierStokesSolver` and `ADRSolver` — not `CompressibleFlowSolver`. This confirms the thesis contribution framing.
- The SVV operator is a local modal-space filter on expansion coefficients (not tied to global matrix assembly), which reduces the coupling cost for porting to the compressible solver.
- Thesis introduction (§1.1) drafted by Valen with Claude feedback; chapter structure reorganized (INTRODUCCION → ESTADO_DEL_ARTE → MOTIVACION split → CAMPANA_NUMERICA merged).
- A 22-entry `.bib` seed file exists with entries marked `[OK]` or `[CHECK]`.
- The project file system takes precedence over memory when they conflict.

**On the horizon**

- **Open SVV design decisions** requiring explicit thesis justification:
  - Normalization scale for SVV coefficient (no physical viscosity in Euler; `h/(P+1)` × characteristic velocity proposed)
  - Choice of variables (conservative, primitive, or entropic)
  - Kernel selection
- **κ(P) measurement campaign**: Protocol established using `IsentropicVortex16` test cases at P=1,3,8; timing via step-difference method; single-core pinned with `taskset`; cross-calibration run needed on ARCHER2 for communication overhead
- **NekMesh/mesh workflow**: NACA 0012 and "Aconcagua" (rocket geometry, highest-priority case) meshes pending. STL geometry limitation identified — faceted, not parametric — with three ranked routes: rebuild as STEP/CAD (preferred), Gmsh linear + spherigon curving, or `mcf` on STL (requires verification of STL CAD backend in 5.10.0)
- **Strong-scaling study** on ARCHER2 separate from local κ measurement
- **Entropy-theoretic framework**: Valen plans to engage with this after consolidating spectral/hp numerical structure understanding — needed for theoretical support but not direct implementation

**Key learnings & principles**

- SVV operates exclusively on high polynomial modes — this directly motivates its design and connects to Gibbs phenomenon mechanics (Dirichlet kernel side lobes, truncation removing modes needed for destructive cancellation)
- Gibbs relative overshoot (~8.95%, Wilbraham constant) does not decay with h- or p-refinement; only the affected region contracts. hp-refinement with mesh aligned to discontinuity avoids it entirely
- The (2P+1) CFL penalty dominates cost at high P/Mach; finest mesh level in a refinement sweep consumes ~94% of total budget; doubling resolution in 3D costs ~16×
- κ exponent b distinguishes sum-factorized from dense-matrix operator implementations — factor ~2 difference at P4 vs P8, with direct consequences for TGV budget extrapolation. The exponent b is machine-portable; the prefactor requires one ARCHER2 cross-calibration
- Benchmark scripts must explicitly *unset* the MPI oversubscription flag and use `--bind-to core` to avoid SMT distortion of cost measurements
- For OpenFOAM comparison: use `simulationType laminar` to make the comparison one of discretization schemes, not physical models; TGV at Re=1600 is the verification case
- Nektar++ `CompressibleFlowSolver` has no explicit turbulence closure — operates as DNS formally, implicit LES when under-resolved, with dissipation from Riemann flux, artificial viscosity (Persson–Peraire sensor), and modal filtering
- Do not assert absence of content in references without verification (Claude made and corrected errors about K&S coverage)

**Approach & patterns**

- Valen drafts all theoretical and prose content first; Claude verifies, stress-tests, and handles mechanical outputs (bibliography entries, scripts, tables). Preferred session shape for entropy-theoretic material: Valen explains arguments to Claude rather than Claude summarizing sources
- Evidence marking system in use: `[V]/[V-bib]/[F]/[X]`; a split into `[V-C]` (Claude-read) and `[V-A]` (author-verified) has been recommended, with only `[V-A]` sources entering thesis text
- Three project files identified as missing and recommended: `DECISIONES.md` (numbered, dated decisions), `BIB_ESTADO.md` (source verification state), `COSTO.md` (computational modeling and infrastructure)
- **Writing style rules** (strict): impersonal voice, dense prose, no transitional filler, no symmetrical hedging, no reflexive three-item lists, no rhetorical bridge questions, no closing summaries. Em-dashes replaced by commas. Technical terms introduced only when context supports them. Register: "el ingeniero de a pie." LaTeX paragraphs separated by blank lines, never `\newline`
- Prefers concise, technically precise responses; uses `/token-efficiency` directive
- Directly corrects Claude when claims are wrong; responses should distinguish environment problems from code problems and flag non-obvious traps proactively

**Tools & resources**

- **Nektar++ 5.10.0** (source build, WSL2 Ubuntu 22.04), `~/nektar-install`
- **OpenFOAM** (finite volume comparison solver)
- **ARCHER2** (target HPC cluster; toolchain: gcc/11.2.0)
- **CMake**, **OpenMPI**, **HDF5**, **FFTW**, **Scotch**, **Boost**, **NekPy**
- **Gmsh** (mesh generation option)
- **WSL2** on Windows with MSYS2/UCRT environment (MinGW confirmed incompatible with Nektar++)
- Machine: `nemepad`, username: `valne`