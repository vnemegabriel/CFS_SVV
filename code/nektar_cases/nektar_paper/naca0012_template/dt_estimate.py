#!/usr/bin/env python3
"""Recommend an explicit RK4 timestep for this mesh + Mach + Reynolds.

  python3 dt_estimate.py [mesh.xml] [conditions.xml]
  python3 dt_estimate.py mesh.xml conditions.xml --mach 1.2 --re 2000 --p 4

Two limits act at once and must be combined:

  dt_conv = Cc * h / ( (U+c) * (2P+1) )        wave transport
  dt_visc = Cv * h^2 / ( nu * (2P+1)^2 )       explicit diffusion
  1/dt_max = 1/dt_conv + 1/dt_visc

Cc = 1.0 is the standard DG estimate.  Cv = 0.091 was CALIBRATED on this case
(NACA0012, P=3, M=0.8, Re=500): measured stable at dt <= 3e-7, NaN at 5e-7,
so dt_max ~ 4e-7.  Re=5000 then predicts 1.7e-6, and 2e-6 was observed stable
for 2000 steps -- i.e. the model is mildly conservative, which is what we want.

At LOW Reynolds the viscous limit dominates and is easy to miss: at Re=500 it
is ~6x more restrictive than the CFL limit.
"""
import argparse, math, sys
import nekcase

Cc, Cv, SAFETY = 1.0, 0.091, 0.5

ap = argparse.ArgumentParser()
ap.add_argument("mesh", nargs="?", default="naca0012.xml")
ap.add_argument("cond", nargs="?", default="conditions.xml")
ap.add_argument("--mach", type=float, help="override Mach")
ap.add_argument("--re",   type=float, help="override Reynolds")
ap.add_argument("--p",    type=int,   help="override polynomial order P (=NUMMODES-1)")
ap.add_argument("--conv-times", type=float, default=27.0, help="target run length")
a = ap.parse_args()

par = nekcase.parameters(a.cond)
ref = nekcase.reference(par)
msh = nekcase.mesh_metrics(a.mesh)
P = a.p if a.p else nekcase.nummodes(a.cond) - 1

# apply overrides consistently: U follows Mach, nu follows Reynolds
M  = a.mach if a.mach else ref["M"]
Re = a.re   if a.re   else ref["Re"]
U  = M * ref["c"]
nu = U * ref["L"] / Re
mu = nu * ref["rho"]

h, np1 = msh["hmin"], 2 * P + 1
lam = U + ref["c"]

dt_conv = Cc * h / (lam * np1)
dt_visc = Cv * h * h / (nu * np1 * np1)
dt_max  = 1.0 / (1.0 / dt_conv + 1.0 / dt_visc)
dt      = SAFETY * dt_max

t_conv  = ref["L"] / U
nsteps  = a.conv_times * t_conv / dt

# worst case artificial viscosity, if ShockCaptureType = Physical
mu0 = par.get("mu0", 1.0)
nu_av = mu0 * (h / max(P, 1)) * lam

print("mesh    : %d elements, h_min = %.4e, h_max = %.4e" % (msh["nelmt"], h, msh["hmax"]))
print("state   : M = %.3f  Re = %.0f  P = %d   U = %.2f m/s  c = %.2f m/s" % (M, Re, P, U, ref["c"]))
print("          mu = %.6g Pa.s   nu = %.6g m2/s" % (mu, nu))
print()
print("dt_conv = %.3e   (CFL / wave transport)" % dt_conv)
print("dt_visc = %.3e   (explicit diffusion)%s" % (dt_visc, "   <-- LIMITING" if dt_visc < dt_conv else ""))
print("dt_max  = %.3e   combined" % dt_max)
print()
print(">>> TimeStep = %.1e            (%.0f%% safety)" % (dt, SAFETY * 100))
print(">>> NumSteps = %d              (%.0f convective times, t_conv = %.3e s)"
      % (round(nsteps, -3), a.conv_times, t_conv))
print()
if nu_av > 0.5 * nu:
    print("WARNING: with ShockCaptureType=Physical the artificial viscosity can reach")
    print("         nu_av ~ %.3g m2/s (mu0=%g), i.e. %.1fx the physical nu." % (nu_av, mu0, nu_av / nu))
    print("         It only fires where the sensor detects compression, so it does not")
    print("         bite until the shock forms -- but if the run dies late, that is why.")
    print("         Fix: lower mu0, or drop TimeStep further, or run with AV Off.")
