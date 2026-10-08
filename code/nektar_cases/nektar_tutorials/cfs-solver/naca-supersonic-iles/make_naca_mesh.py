#!/usr/bin/env python3
"""
Block-structured C-grid around a NACA 4-digit section, extruded into a
spanwise-periodic slab, for supersonic implicit-LES with Nektar++.

    pip install gmsh
    python3 make_naca_mesh.py --level coarse
    python3 make_naca_mesh.py --level les

Output: naca_sup_3d.msh  (Gmsh 2.2 ASCII, curved hexahedra)

------------------------------------------------------------------------
TOPOLOGY
------------------------------------------------------------------------
Four transfinite blocks, the classic C-grid:

        arc_top ______________ out_top          y = +H
          |  wake-upper block   |
          |                     |
   nose  /|--------- TE --------|  out_mid      y = 0
    -R  \ |  upper/lower blocks |
          |_____________________|
        arc_bot                out_bot          y = -H

The outer arc is a circle of radius R centred on the leading edge,
truncated at x = 1 so that H = sqrt(R^2 - 1).  Truncating there (rather
than at the 90-degree point, x = 0) makes the trailing-edge radial line
exactly vertical and the wake-block outer boundary exactly horizontal,
which removes essentially all cell skew from the trailing-edge region --
the one place in a C-grid where skew actually hurts.

------------------------------------------------------------------------
DOMAIN SIZING RULE (why R and LX are what they are)
------------------------------------------------------------------------
At M_inf > 1 every disturbance is confined to a Mach cone of half-angle

    mu = arcsin(1/M)                     (41.81 deg at M = 1.5)

measured from the free-stream direction.  The leading-edge shock is a
little steeper than mu -- it is a weak oblique shock, not a Mach wave --
but mu is the right first-order estimate.  Two consequences:

 1. Nothing propagates upstream outside the cone, so the inflow boundary
    may sit close to the body without polluting the solution.  The C-grid
    puts it at R = 8c essentially for free.

 2. The leading-edge shock reaches height y = H at a streamwise distance
    x = H / tan(mu) downstream of the leading edge.  We want the shock to
    leave through the OUTLET, where PressureOutflow extrapolates cleanly
    once the normal Mach number exceeds one, rather than through the
    far-field boundary, where a Dirichlet free-stream state reflects it
    straight back onto the aerofoil.  That requires

        H / tan(mu) > LX        i.e.        H > LX * tan(mu)

    With H = 7.94c, LX = 8c, M = 1.5:  7.94 > 8 * 0.894 = 7.15.  Satisfied.

    Raising the Mach number shrinks mu, lays the shock flatter, and makes
    the condition easier.  Dropping towards M = 1.1 gives mu -> 65 deg and
    tan(mu) -> 2.14, so you would need H > 2.14 * LX.  The check below is
    enforced at run time -- it will refuse to build a mesh that violates it.

------------------------------------------------------------------------
WALL-NORMAL SPACING
------------------------------------------------------------------------
Y_FIRST is the height of the first ELEMENT, not of the first solution
point.  Nektar++ places NUMMODES+1 Gauss-Lobatto-Legendre points in each
element, and GLL points cluster towards the element ends, so the first
solution point off the wall sits at a fixed fraction of the element:

    P = 2  (4 GLL pts)   y_1 = 0.276 * Y_FIRST
    P = 3  (5 GLL pts)   y_1 = 0.173 * Y_FIRST
    P = 4  (6 GLL pts)   y_1 = 0.117 * Y_FIRST

Raising the polynomial order therefore buys near-wall resolution for
free, without touching the mesh.  Size Y_FIRST from the target y+ using
these factors; README.md does the arithmetic for the chosen Reynolds
number.
"""

import argparse
import math

import gmsh

# --------------------------------------------------------------------------
# Presets.  'coarse' validates the pipeline end to end in minutes.  'les' is
# a genuine wall-resolved slab -- read README.md before launching it.
# --------------------------------------------------------------------------
LEVELS = {
    "coarse": dict(N_AF=60, N_WAKE=25, N_RAD=25, NZ=8, Y_FIRST=2.0e-3),
    "medium": dict(N_AF=120, N_WAKE=45, N_RAD=40, NZ=16, Y_FIRST=5.0e-4),
    "les": dict(N_AF=220, N_WAKE=80, N_RAD=65, NZ=32, Y_FIRST=6.0e-5),
}

# Physical group tags.  NekMesh maps these straight onto Nektar++ composite
# IDs, so C[1] is the fluid volume, C[2] the aerofoil, and so on.  The
# session file depends on these numbers -- do not renumber casually.
TAG_VOLUME = 1
TAG_WALL = 2
TAG_FARFIELD = 3
TAG_OUTFLOW = 4
TAG_PER_LO = 5  # z = 0
TAG_PER_HI = 6  # z = span


def naca4_thickness(x, t):
    """Half-thickness of a NACA 4-digit section.

    The -0.1036 quartic coefficient is the closed-trailing-edge variant.
    The textbook open-TE value (-0.1015) leaves a finite-thickness base,
    which forces a blunt-base separation and demands a far finer wake --
    not what you want here.
    """
    return 5.0 * t * (
        0.2969 * math.sqrt(x)
        - 0.1260 * x
        - 0.3516 * x * x
        + 0.2843 * x**3
        - 0.1036 * x**4
    )


def naca4_camber(x, m, p):
    """Camber line ordinate and slope."""
    if m == 0.0 or p == 0.0:
        return 0.0, 0.0
    if x < p:
        yc = m / (p * p) * (2.0 * p * x - x * x)
        dyc = 2.0 * m / (p * p) * (p - x)
    else:
        yc = m / ((1.0 - p) ** 2) * ((1.0 - 2.0 * p) + 2.0 * p * x - x * x)
        dyc = 2.0 * m / ((1.0 - p) ** 2) * (p - x)
    return yc, dyc


def naca4_surface(digits, n_pts):
    """Cosine-spaced upper and lower surfaces, LE at (0,0), TE at (1,0)."""
    m = int(digits[0]) / 100.0
    p = int(digits[1]) / 10.0
    t = int(digits[2:]) / 100.0

    upper, lower = [], []
    for i in range(n_pts):
        beta = math.pi * i / (n_pts - 1)
        x = 0.5 * (1.0 - math.cos(beta))
        yt = naca4_thickness(x, t)
        yc, dyc = naca4_camber(x, m, p)
        th = math.atan(dyc)
        upper.append((x - yt * math.sin(th), yc + yt * math.cos(th)))
        lower.append((x + yt * math.sin(th), yc - yt * math.cos(th)))
    upper[0] = lower[0] = (0.0, 0.0)
    upper[-1] = lower[-1] = (1.0, 0.0)
    return upper, lower


def radial_growth(y_first, height, n):
    """Geometric ratio marching from y_first to `height` in n elements."""
    target = height / y_first
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(300):
        mid = 0.5 * (lo + hi)
        s = n if abs(mid - 1.0) < 1e-12 else (mid**n - 1.0) / (mid - 1.0)
        if s < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def build(args):
    cfg = LEVELS[args.level]
    n_af, n_wake, n_rad = cfg["N_AF"], cfg["N_WAKE"], cfg["N_RAD"]
    nz, y_first = cfg["NZ"], cfg["Y_FIRST"]

    R, LX, SPAN = args.radius, args.downstream, args.span
    if R <= 1.05:
        raise SystemExit("--radius must exceed the chord.")
    H = math.sqrt(R * R - 1.0)

    mu = math.asin(1.0 / args.mach)
    required = LX * math.tan(mu)
    if H <= required:
        raise SystemExit(
            f"Domain sizing violated.  At M = {args.mach} the Mach angle is "
            f"{math.degrees(mu):.2f} deg, so the leading-edge shock needs a "
            f"far-field height H > LX*tan(mu) = {required:.2f}c to exit "
            f"through the outlet.  This mesh has H = {H:.2f}c.  Increase "
            f"--radius or reduce --downstream."
        )

    growth = radial_growth(y_first, H, n_rad)
    if growth > 1.30:
        print(
            f"  WARNING: radial growth ratio {growth:.3f} > 1.30.  Stretching "
            f"this aggressive degrades the DG diffusion operator and pollutes "
            f"the artificial-viscosity sensor.  Raise N_RAD or Y_FIRST."
        )

    gmsh.initialize()
    gmsh.model.add("naca_supersonic")
    geo = gmsh.model.geo

    upper, lower = naca4_surface(args.naca, 200)

    # ---- aerofoil ---------------------------------------------------------
    p_le = geo.addPoint(0.0, 0.0, 0.0)
    p_te = geo.addPoint(1.0, 0.0, 0.0)
    c_up = geo.addSpline(
        [p_le] + [geo.addPoint(x, y, 0.0) for x, y in upper[1:-1]] + [p_te]
    )
    c_lo = geo.addSpline(
        [p_le] + [geo.addPoint(x, y, 0.0) for x, y in lower[1:-1]] + [p_te]
    )

    # ---- outer boundary ---------------------------------------------------
    p_nose = geo.addPoint(-R, 0.0, 0.0)
    p_arc_top = geo.addPoint(1.0, H, 0.0)
    p_arc_bot = geo.addPoint(1.0, -H, 0.0)
    p_out_top = geo.addPoint(LX, H, 0.0)
    p_out_bot = geo.addPoint(LX, -H, 0.0)
    p_out_mid = geo.addPoint(LX, 0.0, 0.0)

    # Arcs sweep 97.2 deg at R = 8c -- comfortably under gmsh's 180 deg limit.
    a_top = geo.addCircleArc(p_nose, p_le, p_arc_top)
    a_bot = geo.addCircleArc(p_nose, p_le, p_arc_bot)
    l_out_top = geo.addLine(p_arc_top, p_out_top)
    l_out_bot = geo.addLine(p_arc_bot, p_out_bot)
    l_exit_top = geo.addLine(p_out_mid, p_out_top)
    l_exit_bot = geo.addLine(p_out_mid, p_out_bot)

    # ---- internal radial / wake lines -------------------------------------
    l_rad_nose = geo.addLine(p_le, p_nose)
    l_rad_te_top = geo.addLine(p_te, p_arc_top)  # vertical by construction
    l_rad_te_bot = geo.addLine(p_te, p_arc_bot)
    l_wake = geo.addLine(p_te, p_out_mid)

    # ---- four four-sided transfinite blocks -------------------------------
    s_up = geo.addPlaneSurface(
        [geo.addCurveLoop([c_up, l_rad_te_top, -a_top, -l_rad_nose])]
    )
    s_lo = geo.addPlaneSurface(
        [geo.addCurveLoop([c_lo, l_rad_te_bot, -a_bot, -l_rad_nose])]
    )
    s_wk_top = geo.addPlaneSurface(
        [geo.addCurveLoop([l_wake, l_exit_top, -l_out_top, -l_rad_te_top])]
    )
    s_wk_bot = geo.addPlaneSurface(
        [geo.addCurveLoop([l_wake, l_exit_bot, -l_out_bot, -l_rad_te_bot])]
    )
    geo.synchronize()

    # ---- transfinite distributions ----------------------------------------
    # Bump < 1 clusters towards BOTH ends of a curve: at the leading edge
    # (stagnation, strong favourable gradient) and at the trailing edge
    # (wake formation, shock foot).  This is the single most consequential
    # streamwise clustering choice in the mesh.
    for c in (c_up, c_lo, a_top, a_bot):
        geo.mesh.setTransfiniteCurve(c, n_af + 1, "Bump", 0.08)
    for c in (l_rad_nose, l_rad_te_top, l_rad_te_bot, l_exit_top, l_exit_bot):
        geo.mesh.setTransfiniteCurve(c, n_rad + 1, "Progression", growth)
    for c in (l_wake, l_out_top, l_out_bot):
        geo.mesh.setTransfiniteCurve(c, n_wake + 1, "Progression", 1.045)

    geo.mesh.setTransfiniteSurface(s_up, "Left", [p_le, p_te, p_arc_top, p_nose])
    geo.mesh.setTransfiniteSurface(s_lo, "Left", [p_le, p_te, p_arc_bot, p_nose])
    geo.mesh.setTransfiniteSurface(
        s_wk_top, "Left", [p_te, p_out_mid, p_out_top, p_arc_top]
    )
    geo.mesh.setTransfiniteSurface(
        s_wk_bot, "Left", [p_te, p_out_mid, p_out_bot, p_arc_bot]
    )
    for s in (s_up, s_lo, s_wk_top, s_wk_bot):
        geo.mesh.setRecombine(2, s)
    geo.synchronize()

    # ---- spanwise extrusion into a periodic slab --------------------------
    # recombine=True turns the extruded quad layers into hexahedra.  The two
    # z-faces receive identical node distributions by construction, which is
    # exactly what a Nektar++ periodic boundary pair requires.
    ext = []
    for s in (s_up, s_lo, s_wk_top, s_wk_bot):
        ext.append(
            geo.extrude([(2, s)], 0.0, 0.0, SPAN, numElements=[nz], recombine=True)
        )
    geo.synchronize()
    volumes = [t for e in ext for (d, t) in e if d == 3]

    # ---- classify boundary surfaces by geometry ---------------------------
    # Trusting the ordering of extrude()'s return value is brittle; bounding
    # boxes are not.
    bnd = gmsh.model.getBoundary(
        [(3, v) for v in volumes], combined=True, oriented=False
    )
    eps = 1e-6
    wall, farfield, outflow, per_lo, per_hi = [], [], [], [], []
    for dim, tag in bnd:
        if dim != 2:
            continue
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
        if abs(z1 - z0) < eps:
            (per_lo if abs(z0) < eps else per_hi).append(tag)
        elif abs(x0 - LX) < 1e-6 and abs(x1 - LX) < 1e-6:
            outflow.append(tag)
        elif -0.05 < x0 and x1 < 1.05 and max(abs(y0), abs(y1)) < 0.25:
            wall.append(tag)
        else:
            farfield.append(tag)

    print(
        f"  surfaces  wall={len(wall)}  farfield={len(farfield)}  "
        f"outflow={len(outflow)}  per_lo={len(per_lo)}  per_hi={len(per_hi)}"
    )
    if len(wall) != 2 or len(per_lo) != 4 or len(per_lo) != len(per_hi):
        print(
            "  WARNING: unexpected surface counts.  Expected wall=2 (upper and "
            "lower aerofoil), per_lo=per_hi=4 (one per block).  Inspect the "
            "mesh in the gmsh GUI before trusting the boundary conditions."
        )

    gmsh.model.addPhysicalGroup(3, volumes, TAG_VOLUME)
    gmsh.model.setPhysicalName(3, TAG_VOLUME, "fluid")
    for tag, surfs, name in (
        (TAG_WALL, wall, "wall"),
        (TAG_FARFIELD, farfield, "farfield"),
        (TAG_OUTFLOW, outflow, "outflow"),
        (TAG_PER_LO, per_lo, "periodic_lo"),
        (TAG_PER_HI, per_hi, "periodic_hi"),
    ):
        gmsh.model.addPhysicalGroup(2, surfs, tag)
        gmsh.model.setPhysicalName(2, tag, name)

    # ---- curved, high-order mesh ------------------------------------------
    # The aerofoil must be represented to at least the order of the solution
    # expansion.  If it is not, the geometric error dominates and the
    # high-order convergence you are paying for never appears.
    gmsh.option.setNumber("Mesh.ElementOrder", args.order)
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)  # NekMesh reads 2.2

    gmsh.model.mesh.generate(3)
    gmsh.write(args.output)

    nelm = sum(len(gmsh.model.mesh.getElements(3, v)[1][0]) for v in volumes)
    print(f"\n  wrote {args.output}")
    print(f"  hexahedra       : {nelm:,}")
    print(f"  far-field height: {H:.3f} c   (needs > {required:.3f} c)")
    print(f"  Mach angle      : {math.degrees(mu):.2f} deg")
    print(f"  radial growth   : {growth:.4f}")
    print(f"  first element   : {y_first:.2e} c")
    for p in (2, 3, 4):
        print(f"  DOF at P={p}      : {nelm * (p + 1) ** 3 * 5:,}")
    gmsh.finalize()


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--naca", default="0012", help="4-digit designation")
    ap.add_argument("--level", default="coarse", choices=sorted(LEVELS))
    ap.add_argument("--mach", type=float, default=1.5, help="for the sizing check")
    ap.add_argument("--radius", type=float, default=8.0, help="far-field radius / c")
    ap.add_argument("--downstream", type=float, default=8.0, help="outlet x / c")
    ap.add_argument("--span", type=float, default=0.1, help="spanwise extent / c")
    ap.add_argument("--order", type=int, default=4, help="geometric element order")
    ap.add_argument("--output", default="naca_sup_3d.msh")
    build(ap.parse_args())
