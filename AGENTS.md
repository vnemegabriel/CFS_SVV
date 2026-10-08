# TESIS — SVV in Nektar++ CompressibleFlowSolver

Shared instructions for Claude Code and opencode. Harness-specific notes live in `CLAUDE.md` and `opencode.json`.

## Goal

Implement and calibrate a Spectral Vanishing Viscosity (SVV) operator in `CompressibleFlowSolver`, following `code/SVV_nektar_implementation_plan.md` (file map) and `casos/CASOS.md` (validation levels). The filesystem wins over memory when they disagree.

## Two Nektar++ installs (never mix)

| Role | Source | Build | Binaries |
|---|---|---|---|
| **base** (reference, read-only) | `~/nektar` @ `master` | `~/nektar/build` | `/opt/nektar++` |
| **dev** (construction) | `~/nektar` @ `feature/svv-cfs` | `~/nektar/build-svv` (RelWithDebInfo) | `~/nektar/build-svv/dist` |

- Run solvers only via `harness/bin/nk base|dev <cmd>`, never with bare `PATH`.
- Build only via `harness/bin/svv-build`; test via `harness/bin/svv-test quick|full`.
- Forbidden: `make install` in `~/nektar/build`, writing to `/opt/nektar++`, `git push`, editing `~/nektar` outside `feature/svv-cfs`.
- Hooks enforce this. If blocked, change approach; do not retry.

## Where code changes go

Only `~/nektar/solvers/CompressibleFlowSolver/` (incl. `Tests/` and its `CMakeLists.txt`). `library/` is read-only unless a decision is logged in `DECISIONES.md`.

Patterns: `ArtificialDiffusion/NonSmoothShockCapture.{h,cpp}` (factory plugin) and `library/StdRegions/Std*Exp.cpp::v_SVVLaplacianFilter` (kernel). Also cover the `DiffuseCoeffs` path (implicit/ALE).

## Code conventions

- Nektar++ style: repo clang-format, `Array<OneD,...>`, `NekFactory`, session keys via `m_session->LoadParameter`.
- No calibrated constants in C++: kernel, `SVVCutoffRatio`, `SVVDiffCoeff` come from the session XML. A code default needs a justification in `DECISIONES.md`.
- Nothing case-specific (case names, Mach, geometry) inside the solver.
- Small commits on `feature/svv-cfs`, one change each; commit messages in Spanish, imperative.

## Validation (in order; each level gates the next)

1. **Builds**: `svv-build`.
2. **Gate**: `svv-test quick` → own SVV tests + 4 base CFS tests. Before merge: `svv-test full` (all non-parallel CFS).
3. **Transparency**: SVV off ⇒ dev reproduces base L2 norms (`harness/bin/validar.py transparencia`).
4. **Calibration**: `validar.py ronda` on the *evolve* set (`harness/validacion/conjuntos.json`).
5. **Holdout**: `validar.py holdout`. Orchestrator only, at milestones; its numbers are never shown to proposers.

## Evidence marks

`[V-A]` author-verified, `[V-C]` read by Claude, `[V-bib]` bibliographic record only, `[F]` unverified, `[X]` error found. Only `[V-A]` enters thesis text.

## Thesis prose (written in Spanish)

Impersonal voice, dense prose, no transitional filler, no reflexive three-item lists, no rhetorical questions, no closing summaries, commas instead of em-dashes. LaTeX paragraphs separated by a blank line. Valentín writes theory; agents verify and do mechanical work.

## Repo map

- `code/` — implementation plan, 1D prototypes (`cfs1d/`, `svv_spectral_hp_1d/`), loose cases.
- `casos/casos/` — campaign runner (`casos.py`); machine `nemepad` = base, `nemepad-dev` = dev (`CAMP_MAQUINA=nemepad-dev`).
- `memory/` — verified technical notes (shock→SVV chain, CFS call sequence, K&S §10.1).
- `raw/` — sources (papers, books). Read-only.
- `harness/` — scripts, hooks, validation, ledger.
