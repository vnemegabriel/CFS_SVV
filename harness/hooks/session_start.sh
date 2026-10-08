#!/usr/bin/env bash
# SessionStart: estado en ≤4 líneas.
H="$(cd "$(dirname "$0")/.." && pwd)"; NEK="$HOME/nektar"
br=$(git -C "$NEK" branch --show-current 2>/dev/null)
dev="sin build"; [ -x "$NEK/build-svv/dist/bin/CompressibleFlowSolver" ] && dev="build OK"
ch=$(git -C "$NEK" status --porcelain 2>/dev/null | wc -l)
echo "[harness] ~/nektar: rama $br, $ch archivos modificados, dev: $dev"
[ -s "$H/validacion/ledger.jsonl" ] && "$H/bin/validar.py" historial -n 3 2>/dev/null | sed 's/^/[ledger] /'
exit 0
