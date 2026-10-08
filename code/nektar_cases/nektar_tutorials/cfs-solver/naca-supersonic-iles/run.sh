#!/usr/bin/env bash
#
# Supersonic NACA implicit-LES pipeline.
#
#   ./run.sh coarse 4          # level, number of MPI ranks
#   ./run.sh les    256
#
# Stops at the first failure rather than marching on with a broken mesh.

set -euo pipefail

LEVEL="${1:-coarse}"
NPROC="${2:-4}"
CASE="naca_sup"

# --------------------------------------------------------------------------
# 0. Preflight.  This build was configured with NEKTAR_USE_MPI=OFF, which
#    means `mpirun -np N CompressibleFlowSolver` launches N INDEPENDENT
#    SERIAL solvers, each one solving the whole problem and each one writing
#    over the others' output.  It is not a parallel run.  See README.md.
# --------------------------------------------------------------------------
if ! CompressibleFlowSolver --version 2>&1 | grep -qi 'MPI'; then
    echo "WARNING: this CompressibleFlowSolver does not report MPI support."
    echo "         Check 'CompressibleFlowSolver --version'.  If MPI is"
    echo "         absent, rebuild before attempting a 3D run - see README.md."
    echo
fi

# --------------------------------------------------------------------------
# 1. Mesh
# --------------------------------------------------------------------------
echo "==> generating mesh (level=${LEVEL})"
python3 make_naca_mesh.py --level "${LEVEL}" --output "${CASE}_3d.msh"

# --------------------------------------------------------------------------
# 2. Convert to Nektar++ geometry.
#    :xml:uncompress keeps the geometry human-readable so you can inspect
#    the <COMPOSITE> block and confirm the physical-group -> composite map
#    that session_naca_sup.xml assumes.  Drop it once you trust the mesh;
#    compressed geometry is much smaller and loads faster.
# --------------------------------------------------------------------------
echo "==> converting mesh"
NekMesh "${CASE}_3d.msh" "${CASE}_3d.xml:xml:uncompress"

echo
echo "    composites found:"
grep -o '<C ID="[0-9]*"[^>]*>[^<]*' "${CASE}_3d.xml" | head -20 || true
echo
echo "    Expected: C[1] volume, C[2] wall, C[3] farfield,"
echo "              C[4] outflow, C[5] z=0, C[6] z=span."
echo "    If these differ, fix BOUNDARYREGIONS in session_naca_sup.xml."
echo

# --------------------------------------------------------------------------
# 3. Solve
# --------------------------------------------------------------------------
echo "==> solving on ${NPROC} rank(s)"
if [ "${NPROC}" -gt 1 ]; then
    mpirun -np "${NPROC}" CompressibleFlowSolver "${CASE}_3d.xml" session_naca_sup.xml
else
    CompressibleFlowSolver "${CASE}_3d.xml" session_naca_sup.xml
fi

# --------------------------------------------------------------------------
# 4. Post-process.
#    The extra field modules are what make the result readable:
#      shockcapture  - the artificial viscosity field.  ALWAYS look at this.
#                      If it is non-zero anywhere except at the shock, the
#                      LES is being damped and the result is not trustworthy.
#      vorticity     - resolved turbulent structures.
#      gradient      - for Q-criterion / lambda2 isosurfaces.
# --------------------------------------------------------------------------
echo "==> post-processing"
FieldConvert -m vorticity \
    "${CASE}_3d.xml" session_naca_sup.xml "${CASE}.fld" "${CASE}_flow.vtu"

if [ -d "${CASE}_avg.fld" ] || [ -f "${CASE}_avg.fld" ]; then
    FieldConvert "${CASE}_3d.xml" session_naca_sup.xml \
        "${CASE}_avg.fld" "${CASE}_mean.vtu"
fi

echo
echo "done.  open ${CASE}_flow.vtu in ParaView."
echo "forces are in ${CASE}_forces.fce"
