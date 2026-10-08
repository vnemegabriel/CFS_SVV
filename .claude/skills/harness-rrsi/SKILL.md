---
name: harness-rrsi
description: Regularized self-improvement of the harness itself (AGENTS.md, agent prompts, hooks, scripts) driven by the ledger and observed failures. Use when an agent fails repeatedly, token cost grows, or the user asks to revise the harness.
---
# Harness self-improvement (RRSI on the harness)

Evolving object: `AGENTS.md`, `CLAUDE.md`, `.claude/agents/*.md`, `harness/opencode/*.md`, `harness/hooks/*`, `harness/bin/*`. Models are fixed.

1. **Evidence.** Collect `validar.py historial -n 20`, `harness/logs/` (build/test failures, `oc-*.err`), `harness/validacion/harness_ledger.jsonl`. Count failures per agent, `oc` retries escalated to Haiku/Sonnet, budget/leak rejections.
2. **Proposal.** One edit per iteration (budget fixed at 1) with a falsifiable hypothesis ("the coder skips DiffuseCoeffs because the prompt never names it"). Nothing specific to a case or a round.
3. **Critique.** `oc svv-critic` on the harness diff: reject rules that memorize a case, steps nobody runs, text that only adds tokens.
4. **Acceptance.** Compare the next 3 equivalent tasks with the previous 3: accept if failures or tokens drop without the other rising. Within noise, prefer the shorter version.
5. **Pruning.** Any instruction not cited and not preventing a failure in the last 10 tasks is a deletion candidate.
6. **Log.** One line in `harness/validacion/harness_ledger.jsonl`: `{"date","file","hypothesis","evidence","verdict"}`.

Changes to hooks or permissions always need user confirmation.
