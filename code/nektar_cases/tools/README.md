# Nektar++ postprocessing

Cross-platform (Windows via WSL / native Linux). Requires FieldConvert on the
WSL/Linux PATH, and `pip install pyvista imageio-ffmpeg` for the viewer.

    cd code/nektar_cases/nektar_tutorials/cfs-solver/cfs-Euler-01

    # all naca_*.chk + naca.fld -> vtu/ with a time-series naca.pvd
    python ../../../tools/nek_post.py --mesh naca.xml \
        --session session_naca_ex-Complete.xml --order 4

    # derived compressible fields
    python ../../../tools/nek_post.py --mesh naca.xml \
        --session session_naca_ex-Complete.xml --modules mach pressure vorticity

    # view / animate
    python ../../../tools/nek_view.py vtu/naca.pvd --field u
    python ../../../tools/nek_view.py vtu/naca.pvd --field u --movie wake.mp4

Open `vtu/naca.pvd` directly in ParaView to get the whole run as one animatable
dataset (works identically on Windows and Linux).

Options
- `--order N` : high-order VTU (VTK Lagrange cells). Curved elements stay curved;
  no visual faceting. ParaView >= 5.9.
- `--npts N`  : alternative — linear VTU subdivided at N equispaced points per
  direction. Bigger files, works with older readers. Do not combine with --order.
- `--modules` : any FieldConvert module (mach, pressure, vorticity, QCriterion,
  wss, ...). Applied to every checkpoint.

## If `--order` warns "Unrecognised config option highorder"

The VTU writer option name is build-dependent. List what your binary accepts:

    FieldConvert -l | grep -i -A8 'vtk\|vtu'

If `highorder` is absent, the build has no Lagrange-cell writer (VTK linked
without it). Use `--npts 6` instead — same smoothness on screen, larger files.
