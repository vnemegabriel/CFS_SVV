#!/usr/bin/env bash
# Stop: si ~/nektar cambió desde el último verde, compila y corre svv-test quick.
# Silencioso si no hay cambios o si todo pasa. Detecta también cambios hechos por opencode.
H="$(cd "$(dirname "$0")/.." && pwd)"; NEK="$HOME/nektar"; ST="$H/state/last_green"
in=$(cat)
[ "$(echo "$in" | jq -r '.stop_hook_active // false')" = true ] && exit 0
[ "$(git -C "$NEK" branch --show-current 2>/dev/null)" = feature/svv-cfs ] || exit 0
[ -f "$NEK/build-svv/CMakeCache.txt" ] || exit 0          # dev todavía no configurado
sig=$( { git -C "$NEK" rev-parse HEAD; git -C "$NEK" diff HEAD;
         git -C "$NEK" ls-files -o --exclude-standard -z | xargs -0 -r cat; } 2>/dev/null | sha1sum | cut -c1-16)
[ -f "$ST" ] && [ "$(cat "$ST")" = "$sig" ] && exit 0
if ! out=$("$H/bin/svv-build" 2>&1); then
  { echo "Gate: el build de dev falla. Corregir antes de terminar:"; echo "$out" | tail -25; } >&2; exit 2
fi
if ! out=$("$H/bin/svv-test" quick 2>&1); then
  { echo "Gate: svv-test quick falla. Corregir o explicar al usuario:"; echo "$out" | tail -20; } >&2; exit 2
fi
echo "$sig" > "$ST"; exit 0
