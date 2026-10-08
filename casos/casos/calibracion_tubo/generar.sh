#!/bin/bash
# Copia la malla y las sesiones del tubo de choque 2D a plantilla/ con las marcas @nombre@.
#
#     ./generar.sh [raíz de nektar]          por defecto ~/nektar-master
#
# Las sesiones originales son pruebas de humo (1 a 10 pasos, Euler explícito de orden 1).
# Se marcan pasos, paso de tiempo e integrador. La implícita conserva su Euler hacia atrás.
set -eu
SRC=${1:-$HOME/nektar-master}/solvers/CompressibleFlowSolver/Tests
DST=$(cd "$(dirname "$0")" && pwd)/plantilla
mkdir -p "$DST"
cp "$SRC/ShockTube_2D_mixedMesh.xml" "$DST/"

for f in "$SRC"/ShockTube_2D_mixedMesh_*.xml; do
    s=$(basename "$f" .xml)
    s=${s#ShockTube_2D_mixedMesh_}
    args=(-e 's|<P> *NumSteps *=[^<]*</P>|<P> NumSteps     = @steps@ </P>|'
          -e 's|<P> *TimeStep *=[^<]*</P>|<P> TimeStep     = @dt@ </P>|'
          -e 's|<P> *IO_InfoSteps *=[^<]*</P>|<P> IO_InfoSteps = @info_steps@ </P>|')
    if ! grep -q "Backward" "$f"; then
        args+=(-e 's|<METHOD>[^<]*</METHOD>|<METHOD> @method@ </METHOD>|'
               -e 's|<ORDER>[^<]*</ORDER>|<ORDER> @order@ </ORDER>|'
               -e '/<VARIANT>/d')
    fi
    sed "${args[@]}" "$f" > "$DST/$s.xml"
done

ls "$DST" | wc -l
grep -o "@[a-z_]*@" "$DST"/*.xml | sort | uniq -c | awk '{print $1, $2}' | sed 's|.*/||' | sort | uniq -c
