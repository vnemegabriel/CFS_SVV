@AGENTS.md

## Model routing (Claude Code)

- **Main session = Opus**: plans, decides, reviews diffs, triggers holdout, writes `DECISIONES.md`. Does not write long C++.
- **Sonnet** (`svv-implementer`, `svv-critic`): bounded C++ changes, diff review, critical reading of results.
- **Haiku** (`nektar-explorer`, `bulk-runner`): searches in `~/nektar`, runs, sweeps, log parsing.
- **opencode (free models)** via `harness/bin/oc <agent> "<task>"`: boilerplate, `.tst/.xml` tests, first-draft code, second-opinion critic. Its output is always reviewed before acceptance.

Cost rule: anything mechanical and verifiable by build/test goes to `oc` first; after two failures, escalate to Haiku/Sonnet.

Reply to the user in English (thesis prose stays Spanish).

## Flow

Implementation + calibration loop: skill `/svv-ciclo`. Harness self-improvement: skill `/harness-rrsi`.

The `Stop` hook builds and runs `svv-test quick` whenever `~/nektar` changed since the last green state (including opencode edits). On failure the trimmed error comes back; fix it before finishing.
