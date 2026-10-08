#!/bin/bash
# Convert Nektar++ checkpoints to VTU with primitive variables, for ParaView.
#
#   ./viz.sh                                  viscous case, every checkpoint
#   ./viz.sh 5                                only checkpoint 5
#   COND=conditions-euler.xml ./viz.sh        Euler case
#   KEEP=1 ./viz.sh                           do not delete previous output
#
# Output is tagged by solver ("ns" or "euler"), read from EQType in the
# conditions file.  Each run deletes ITS OWN previous output first, so stale
# checkpoints never linger, but the other solver's results are left alone.
#
# Reference values (gamma, pInf, q_inf) come from the conditions file via
# nekcase.py, so this stays correct at any Mach / Reynolds without editing.
set -e

MESH=${MESH:-naca0012.xml}
COND=${COND:-conditions.xml}
SESS=$(basename "$MESH" .xml)
OUT=${OUT:-vtk}

read -r TAG G PINF Q MACH RE ALPHA EQ <<<"$(python3 -c "
import nekcase
r = nekcase.reference(nekcase.parameters('$COND'))
print(nekcase.tag('$COND'), r['gamma'], r['p'], r['q'],
      r['M'], r['Re'], r['alpha'],
      nekcase.solverinfo('$COND').get('EQType','?'))")"
GM1=$(python3 -c "print($G-1)")

printf 'solver   : %s  (tag "%s")\n' "$EQ" "$TAG"
printf 'state    : M=%.3f  alpha=%g deg' "$MACH" "$ALPHA"
[ "$TAG" = "ns" ] && printf '  Re=%.0f' "$RE"
printf '\nreference: gamma=%s  pInf=%s  q_inf=%s\n' "$G" "$PINF" "$Q"

mkdir -p "$OUT"
if [ -z "$KEEP" ]; then
  old=$(ls "$OUT"/${TAG}_* 2>/dev/null | wc -l)
  if [ "$old" -gt 0 ]; then
    echo "cleaning  : removing $old stale ${TAG}_* file(s) from $OUT/"
    rm -f "$OUT"/${TAG}_*
  fi
fi

KE="0.5*(rhou*rhou+rhov*rhov)/rho"
P="$GM1*(E-$KE)"                       # static pressure from conservative vars
V2="(rhou*rhou+rhov*rhov)/(rho*rho)"

PRIM=(
  -m fieldfromstring:fieldstr="rhou/rho":fieldname="u"
  -m fieldfromstring:fieldstr="rhov/rho":fieldname="v"
  -m fieldfromstring:fieldstr="$P":fieldname="p"
  -m fieldfromstring:fieldstr="($P-$PINF)/$Q":fieldname="Cp"
  -m fieldfromstring:fieldstr="($P)/(287.058*rho)":fieldname="T"
  -m fieldfromstring:fieldstr="sqrt($V2)/sqrt($G*($P)/rho)":fieldname="Mach"
)

if [ -n "$1" ]; then
  list="${SESS}_$1.chk"
else
  list=$(ls -d ${SESS}_[0-9]*.chk 2>/dev/null | grep -v bak | sort -V)
fi
[ -z "$list" ] && { echo "no ${SESS}_N.chk found"; exit 1; }

echo
for c in $list; do
  n=$(echo "$c" | sed "s/${SESS}_//;s/\.chk//")
  nn=$(printf '%04d' "$n")            # zero-padded: keeps the ParaView series ordered
  echo ">>> $c -> ${TAG}_${nn}.vtu"
  # -e = equispaced output; without it P>1 elements render as flat linear cells
  FieldConvert -f -e "${PRIM[@]}" "$MESH" "$COND" "$c" "$OUT/${TAG}_${nn}.vtu" >/dev/null 2>&1
  # boundary 0 = aerofoil.  NOTE: extract drops derived fields, so this file
  # carries only rho/rhou/rhov/E -- compute Cp in ParaView with a Calculator.
  FieldConvert -f -e -m extract:bnd=0 "$MESH" "$COND" "$c" "$OUT/${TAG}_surf_${nn}.vtu" >/dev/null 2>&1
done

# manifest, so it is always clear what a set of files actually contains
{
  echo "solver     : $EQ"
  echo "conditions : $COND"
  echo "mesh       : $MESH"
  echo "Mach       : $MACH"
  echo "alpha      : $ALPHA deg"
  [ "$TAG" = "ns" ] && echo "Reynolds   : $RE"
  echo "pInf       : $PINF"
  echo "q_inf      : $Q"
  echo "generated  : $(date '+%Y-%m-%d %H:%M')"
} > "$OUT/${TAG}_INFO.txt"

echo
echo "done -> $OUT/"
echo "  volume : $OUT/${TAG}_....vtu       (u v p T Cp Mach + conservative)"
echo "  surface: $OUT/${TAG}_surf_...._b0.vtu"
echo "  info   : $OUT/${TAG}_INFO.txt"
echo
echo "  ParaView surface Cp Calculator:"
echo "    ($GM1*(E-0.5*(rhou^2+rhov^2)/rho) - $PINF) / $Q"
