#!/usr/bin/env python3
"""Parametric NACA00xx O-grid for Nektar++, built with the Gmsh Python API.

GRADING -- the two directions are controlled independently:

  radial (wall-normal)    --first-cell H    height of the first cell at the wall
                          --re RE           pick --first-cell to resolve Re
                          --radial-ratio R  set the growth ratio directly instead

  streamwise (along the   --le-cell H       cell length at the leading and
  aerofoil surface)                         trailing edge
                          --af-ratio R      set the growth ratio directly instead

You give a physical size, the script solves for the geometric progression that
hits it. Pass a ratio instead if you would rather drive it that way. Either way
--report prints the resulting distribution, cell by cell, so the grading is
visible before you mesh anything.

  python3 make_mesh.py --report                        # see the default grading
  python3 make_mesh.py --re 5000 --le-cell 1e-3        # size both directions
  python3 make_mesh.py --radial-ratio 1.08 --af-ratio 1.03
  python3 make_mesh.py --n-af 160 --n-rad 60 --order 4 -o fine.msh

Then:
  NekMesh fine.msh fine.xml:xml:uncompress
  # delete the <EXPANSIONS> block NekMesh adds -- conditions.xml owns it

Structure: four transfinite blocks (upper/lower x fore/aft), quad-recombined,
giving a structured O-grid. Physical groups match the case files:
    1 = domain, 5900 = aerofoil, 5901 = farfield
"""
import argparse, math

ap = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("-o", "--out", default="naca0012_gen.msh")
ap.add_argument("--thickness", type=float, default=0.12, help="NACA 4-digit t/c")
ap.add_argument("--chord", type=float, default=1.0)
ap.add_argument("--r-far", type=float, default=15.0, help="farfield radius, chords")
ap.add_argument("--n-af", type=int, default=100, help="points per surface half")
ap.add_argument("--n-rad", type=int, default=48, help="points in the radial direction")
ap.add_argument("--order", type=int, default=3, help="geometric element order")
# --- grading, radial
ap.add_argument("--first-cell", type=float, default=None, help="wall-normal first cell / chord")
ap.add_argument("--re", type=float, default=None, help="target Re -> sets --first-cell")
ap.add_argument("--radial-ratio", type=float, default=None, help="radial growth ratio (overrides)")
# --- grading, streamwise
ap.add_argument("--le-cell", type=float, default=None, help="streamwise cell at LE/TE / chord")
ap.add_argument("--af-ratio", type=float, default=None, help="streamwise growth ratio (overrides)")
# ---
ap.add_argument("--report", action="store_true", help="print the cell distributions")
ap.add_argument("--gui", action="store_true", help="open the Gmsh GUI when done")
a = ap.parse_args()

# ---------------------------------------------------------------- profile
def naca_half(x, t):
    """Upper surface of a symmetric NACA 4-digit section. Closed trailing edge."""
    return 5 * t * (0.2969 * math.sqrt(x) - 0.1260 * x - 0.3516 * x**2
                    + 0.2843 * x**3 - 0.1036 * x**4)

_np = 2 * a.n_af - 1                       # points per surface, split at the midpoint
xs = [0.5 * (1 - math.cos(math.pi * i / (_np - 1))) for i in range(_np)]

c = a.chord
_pts = [(x * c, naca_half(x, a.thickness) * c) for x in xs]
_mid = a.n_af - 1
# arc length of one surface half (LE -> mid-chord), used to size the streamwise grading
L_af = sum(math.dist(_pts[i], _pts[i + 1]) for i in range(_mid))
L_rad = a.r_far - 0.5 * c

# ---------------------------------------------------------------- grading
def solve_ratio(L, h1, n_int):
    """Geometric ratio r with h1*(r^n - 1)/(r - 1) = L. r = 1 if h1 already fits."""
    if h1 * n_int >= L:
        return 1.0
    lo, hi = 1.0 + 1e-12, 5.0
    for _ in range(300):
        r = 0.5 * (lo + hi)
        if h1 * (r**n_int - 1) / (r - 1) < L:
            lo = r
        else:
            hi = r
    return 0.5 * (lo + hi)

def first_cell(L, r, n_int):
    """Inverse: the first cell that a given ratio produces."""
    return L / n_int if abs(r - 1) < 1e-12 else L * (r - 1) / (r**n_int - 1)

n_rad_int, n_af_int = a.n_rad - 1, a.n_af - 1

# radial
if a.radial_ratio:
    r_rad = a.radial_ratio
    h_rad = first_cell(L_rad, r_rad, n_rad_int)
else:
    if a.first_cell is None:
        a.first_cell = (5.0 / math.sqrt(a.re)) / 8.0 if a.re else 5.0e-3
    h_rad = a.first_cell * c
    r_rad = solve_ratio(L_rad, h_rad, n_rad_int)

# streamwise
if a.af_ratio:
    r_af = a.af_ratio
    h_af = first_cell(L_af, r_af, n_af_int)
else:
    if a.le_cell is None:
        a.le_cell = L_af / n_af_int / c          # uniform by default
    h_af = a.le_cell * c
    r_af = solve_ratio(L_af, h_af, n_af_int)

def cells(L, h1, r, n):
    return [h1 * r**i for i in range(n)] if abs(r - 1) > 1e-12 else [L / n] * n

cr, ca = cells(L_rad, h_rad, r_rad, n_rad_int), cells(L_af, h_af, r_af, n_af_int)

print("NACA00%02d   chord=%g   farfield=%g c   order %d" %
      (round(a.thickness * 100), c, a.r_far, a.order))
print("grid       %d x %d per block  ->  %d elements"
      % (a.n_af, a.n_rad, 4 * n_af_int * n_rad_int))
print()
print("  direction    n     ratio     first        last       total")
print("  radial     %4d   %7.5f   %.4e  %.4e  %.4f" % (a.n_rad, r_rad, cr[0], cr[-1], sum(cr)))
print("  streamwise %4d   %7.5f   %.4e  %.4e  %.4f" % (a.n_af, r_af, ca[0], ca[-1], sum(ca)))
if a.re:
    d99 = 5.0 / math.sqrt(a.re)
    n_in = sum(1 for k in range(len(cr)) if sum(cr[:k + 1]) <= d99 * c)
    print("\n  Re = %.0f  ->  delta99 = %.4g c, %d cells inside it" % (a.re, d99, n_in))
if a.report:
    print("\n  radial cell heights (wall -> farfield):")
    for i in (0, 1, 2, 3, 4, n_rad_int // 2, n_rad_int - 2, n_rad_int - 1):
        print("      cell %3d   h = %.5e   y = %.5e" % (i, cr[i], sum(cr[:i + 1])))
    print("  streamwise cell lengths (LE -> mid-chord):")
    for i in (0, 1, 2, 3, 4, n_af_int // 2, n_af_int - 2, n_af_int - 1):
        print("      cell %3d   s = %.5e   arc = %.5e" % (i, ca[i], sum(ca[:i + 1])))
print()

# ---------------------------------------------------------------- geometry
import gmsh
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 0)
gmsh.model.add("naca")
g = gmsh.model.geo

p_le = g.addPoint(0.0, 0.0, 0.0)
p_te = g.addPoint(c, 0.0, 0.0)

def build(sign):
    ids = ([p_le]
           + [g.addPoint(x, sign * y, 0.0) for x, y in _pts[1:-1]]
           + [p_te])
    return g.addSpline(ids[:_mid + 1]), g.addSpline(ids[_mid:]), ids[_mid]

s1, s2, p_mu = build(+1)                           # upper: LE->mid, mid->TE
s3, s4, p_ml = build(-1)                           # lower

cx, cy, R = 0.25 * c, 0.0, a.r_far                 # four 90-deg arcs: the geo
p_ctr = g.addPoint(cx, cy, 0.0)                    # kernel cannot do 180
p_fl = g.addPoint(cx - R, cy, 0.0)
p_ftop = g.addPoint(cx, cy + R, 0.0)
p_ft = g.addPoint(cx + R, cy, 0.0)
p_fbot = g.addPoint(cx, cy - R, 0.0)
A1 = g.addCircleArc(p_fl, p_ctr, p_ftop)
A2 = g.addCircleArc(p_ftop, p_ctr, p_ft)
A3 = g.addCircleArc(p_ft, p_ctr, p_fbot)
A4 = g.addCircleArc(p_fbot, p_ctr, p_fl)

R_le = g.addLine(p_le, p_fl)
R_mu = g.addLine(p_mu, p_ftop)
R_te = g.addLine(p_te, p_ft)
R_ml = g.addLine(p_ml, p_fbot)

blocks = [
    ([s1, R_mu, -A1, -R_le], [p_le, p_mu, p_ftop, p_fl]),
    ([s2, R_te, -A2, -R_mu], [p_mu, p_te, p_ft, p_ftop]),
    ([s3, R_ml,  A4, -R_le], [p_le, p_ml, p_fbot, p_fl]),
    ([s4, R_te,  A3, -R_ml], [p_ml, p_te, p_ft, p_fbot]),
]
surfs = [g.addPlaneSurface([g.addCurveLoop(loop)]) for loop, _ in blocks]
g.synchronize()

# ---------------------------------------------------------------- structure
# Streamwise: curves running LE->mid grow (ratio r), mid->TE shrink (1/r), so
# cells are fine at BOTH the leading and the trailing edge.
# NOTE: setTransfiniteCurve redistributes points along the spline, so the
# progression below -- not the cosine spacing of the control points -- is what
# actually sets the surface distribution.
for l in (s1, s3):
    gmsh.model.mesh.setTransfiniteCurve(l, a.n_af, "Progression", r_af)
for l in (s2, s4):
    gmsh.model.mesh.setTransfiniteCurve(l, a.n_af, "Progression", 1.0 / r_af)
for l in (A1, A4):                                  # match on the farfield side
    gmsh.model.mesh.setTransfiniteCurve(l, a.n_af, "Progression", r_af)
for l in (A2, A3):
    gmsh.model.mesh.setTransfiniteCurve(l, a.n_af, "Progression", 1.0 / r_af)
# Radial: grows away from the wall.
for l in (R_le, R_mu, R_te, R_ml):
    gmsh.model.mesh.setTransfiniteCurve(l, a.n_rad, "Progression", r_rad)

for s, (_, corners) in zip(surfs, blocks):
    gmsh.model.mesh.setTransfiniteSurface(s, "Left", corners)
    gmsh.model.mesh.setRecombine(2, s)

gmsh.model.addPhysicalGroup(2, surfs, 1)
gmsh.model.addPhysicalGroup(1, [s1, s2, s3, s4], 5900)
gmsh.model.addPhysicalGroup(1, [A1, A2, A3, A4], 5901)
gmsh.model.setPhysicalName(2, 1, "domain")
gmsh.model.setPhysicalName(1, 5900, "aerofoil")
gmsh.model.setPhysicalName(1, 5901, "farfield")

# ---------------------------------------------------------------- mesh
gmsh.option.setNumber("Mesh.RecombineAll", 1)
gmsh.option.setNumber("Mesh.ElementOrder", a.order)
gmsh.option.setNumber("Mesh.HighOrderOptimize", 0)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)     # NekMesh reads 2.2
gmsh.model.mesh.generate(2)
gmsh.model.mesh.setOrder(a.order)
gmsh.write(a.out)

print("wrote %s" % a.out)
print("  NekMesh %s %s:xml:uncompress" % (a.out, a.out.replace(".msh", ".xml")))
if a.gui:
    gmsh.fltk.run()
gmsh.finalize()
