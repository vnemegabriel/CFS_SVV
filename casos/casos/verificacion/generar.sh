#!/bin/bash
# Copia las sesiones del test suite de Nektar++ a plantilla/ y les pone las marcas @nombre@.
#
#     ./generar.sh [raíz de nektar]          por defecto ~/nektar-master
#
# V2 (SquareDomain) es EulerCFE en el original: se porta a Navier–Stokes (CASOS.md).
set -eu
SRC=${1:-$HOME/nektar-master}/solvers/CompressibleFlowSolver/Tests
DST=$(cd "$(dirname "$0")" && pwd)/plantilla
mkdir -p "$DST"

for s in SquareDomain_Euler_2D_AxialFlow SquareDomain_Euler_2D_DiagonalFlow; do
    sed -e 's/VALUE="EulerCFE"/VALUE="@eq@"/' \
        -e '/PROPERTY="EQType"/a\            <I PROPERTY="DiffusionType"         VALUE="@difusion@"          />\n            <I PROPERTY="ViscosityType"         VALUE="Constant"            />' \
        -e '/<PARAMETERS>/a\            <P> mu              = @mu@                        </P>\n            <P> Pr              = @Pr@                        </P>' \
        -e 's/NUMMODES="[0-9]*"/NUMMODES="@nummodes@"/' \
        "$SRC/$s.xml" > "$DST/$s.xml"
done

for s in Couette_WeakDG_LDG_SEM MMS_Compressible_Poiseuille_testIP; do
    sed -e 's/NUMMODES="[0-9]*"/NUMMODES="@nummodes@"/' "$SRC/$s.xml" > "$DST/$s.xml"
done

grep -o "@[a-z_]*@" "$DST"/*.xml | sort | uniq -c
