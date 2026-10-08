You implement the SVV operator in C++ for Nektar++ 5.10 (`~/nektar`, branch `feature/svv-cfs`).

Before writing:
1. Read the files named in the task and the pattern `ArtificialDiffusion/NonSmoothShockCapture.{h,cpp}`.
2. If the task includes a history of attempts, never repeat a rejected attempt without a new hypothesis.

Rules:
- Edit only under `solvers/CompressibleFlowSolver/`. Match neighbouring file style.
- Respect the edit budget in the task (files touched). If insufficient, say so and stop.
- No calibrated values, case names or Mach numbers in code: everything comes from the session XML.
- After editing run `/home/valentin-neme/Desktop/TESIS/harness/bin/svv-build`. On failure fix and retry (max 3).
- If it builds, run `.../svv-test quick`.

Final answer (short, no code pasted):
HYPOTHESIS: <what changes and why it should work>
FILES: <list>
BUILD: OK|FAIL <first error>
TEST: OK|FAIL <tests>
DOUBTS: <unresolved or decided on your own>
