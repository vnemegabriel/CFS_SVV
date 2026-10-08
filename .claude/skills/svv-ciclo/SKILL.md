---
name: svv-ciclo
description: Implementation and calibration loop for the SVV filter in CompressibleFlowSolver with regularized (RRSI) validation. Use when advancing a plan step, calibrating SVV parameters, or deciding whether a change is accepted.
---
# SVV loop

The main session (Opus) orchestrates. Each phase names who executes.

## Phase A — implementation (until `"SVV"` exists in the factory)

Plan steps, one per iteration, in order:
1. `SVVDiffusion.{h,cpp}` registered in `ArtificialDiffusionFactory`, key reading, `GetFluxVector` with the modal kernel.
2. `DiffuseCoeffs` path (implicit/ALE).
3. Selection from the session (`CompressibleFlowSystem.cpp`) and `CMakeLists.txt`.
4. Own tests in `Tests/` (scaffold): constant state (SVV of a constant = 0), Couette with SVV ≈ without, Poiseuille MMS. Names contain `SVV` so `svv-test quick` picks them up.

Per iteration:
- Context: `nektar-explorer` (Haiku) when code must be located.
- Code: `svv-implementer` (Sonnet), which first delegates to `oc svv-coder`.
- Review: `svv-critic` (Sonnet + `oc svv-critic`). LEAK/RISK ⇒ back to implementer with reasons.
- Gate: Stop hook builds and runs `svv-test quick`. Before committing, `svv-test full` via `bulk-runner`.
- Commit (orchestrator proposes, user confirms) with the hypothesis in the message.

End of phase A: `validar.py transparencia` must be OK.

## Phase B — calibration (kernel, `SVVCutoffRatio`, `SVVDiffCoeff`, h/p scaling)

Setup (once): `bulk-runner` generates `harness/validacion/sesiones/*_svv.xml` with `oc svv-scaffold` and resolves PENDIENTE values in `conjuntos.json` (shock-tube final time and step). Then `validar.py ruido` fixes incumbent and τ.

Per round:
1. `validar.py presupuesto`, `historial -n 8`, `podar` → orchestrator picks ONE hypothesis within budget; if `podar` reports stagnation, target an unexplored component; components flagged for pruning are proposed for removal.
2. If the hypothesis touches C++: short phase A (implementer → critic → gate).
3. `bulk-runner`: `validar.py ronda --hipotesis "..." --componentes a,b --params K=V ... [--critico]`.
4. ACCEPT becomes the new incumbent. REJECT/LEAK stays in the ledger as evidence; never retried without a new hypothesis.
5. Every 4 accepted rounds, or before any result goes into the thesis: `validar.py holdout --base` (orchestrator only). If holdout worsens while evolve improves, the search is overfitting: stop, log in `DECISIONES.md`, review the evolve set with the user.

## What goes into the thesis

Only holdout results or `casos/` campaign results with the incumbent frozen. Every parameter choice goes to `DECISIONES.md` with date, ledger round and evidence.
