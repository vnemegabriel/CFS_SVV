#!/usr/bin/env python3
"""Nektar++ AeroForces .fce  ->  Cl / Cd  (+ forces.png).

  python3 forces.py [forces.fce] [conditions.xml]

Reference quantities (rho, U, chord, alpha, q_inf) are read from conditions.xml,
so this stays correct at any Mach or Reynolds number without editing.

The AeroForces filter reports forces along the Direction1/Direction2 axes
declared in the filter (here x and y).  Lift and drag are those rotated into
wind axes by the angle of attack.
"""
import math, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nekcase

fce  = sys.argv[1] if len(sys.argv) > 1 else "forces.fce"
cond = sys.argv[2] if len(sys.argv) > 2 else "conditions.xml"
r = nekcase.reference(nekcase.parameters(cond))
q, a = r["q"], math.radians(r["alpha"])

d = np.atleast_2d(np.loadtxt(fce, comments="#"))
if d.size == 0:
    sys.exit("no data rows in %s" % fce)
t, Fx, Fy = d[:, 0], d[:, 3], d[:, 6]            # time, F1-total, F2-total
Cd = ( Fx * math.cos(a) + Fy * math.sin(a)) / q
Cl = (-Fx * math.sin(a) + Fy * math.cos(a)) / q
tc = t * r["U"] / r["L"]                         # convective times

print("M = %.3f  Re = %.0f  alpha = %g deg   q_inf = %.1f Pa   rows = %d"
      % (r["M"], r["Re"], r["alpha"], q, len(t)))
print("final : t = %.4e s (%.2f conv times)   Cl = %+.5f   Cd = %+.5f"
      % (t[-1], tc[-1], Cl[-1], Cd[-1]))
if len(t) > 10:
    h = slice(len(t) // 2, None)
    print("mean over 2nd half :  Cl = %+.5f +/- %.5f   Cd = %+.5f +/- %.5f"
          % (Cl[h].mean(), Cl[h].std(), Cd[h].mean(), Cd[h].std()))
    if tc[-1] < 10:
        print("NOTE: only %.1f convective times -- still in the start transient." % tc[-1])

fig, ax = plt.subplots(2, 1, figsize=(7, 5), sharex=True)
ax[0].plot(tc, Cl);             ax[0].set_ylabel("Cl"); ax[0].grid(alpha=.3)
ax[1].plot(tc, Cd, color="C1"); ax[1].set_ylabel("Cd"); ax[1].grid(alpha=.3)
ax[1].set_xlabel("convective times   t*U/c")
ax[0].set_title("NACA0012   M=%.2f   Re=%.0f   alpha=%g deg" % (r["M"], r["Re"], r["alpha"]))
fig.tight_layout(); fig.savefig("forces.png", dpi=130)
print("wrote forces.png")
