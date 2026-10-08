// =============================================================================
// cns1d_svv  --  1D compressible Euler / Navier-Stokes with SVV (C0 Galerkin)
//
// Cases
//   --case entropywave   smooth density wave, periodic, exact solution known.
//                        Use it to check spectral accuracy with and without SVV.
//   --case sod           Sod shock tube on [0,1], T = 0.2 (exact Riemann
//                        solution in scripts/exact_riemann.py).
//   --case shuosher      Shu-Osher shock/entropy-wave interaction, [-5,5], T=1.8.
//
// Typical runs
//   ./cns1d_svv --case entropywave --nel 8 --P 8
//   ./cns1d_svv --case sod --nel 40 --P 8 --mu 1e-3                    # NS, no SVV
//   ./cns1d_svv --case sod --nel 40 --P 8 --mu 1e-3 --svv --eps 0.05 --scaled
//   ./cns1d_svv --case sod --nel 40 --P 8 --svv --eps 0.05 --scaled   # Euler + SVV only
//
// SVV amplitude: "--eps e" is a constant eps unless "--scaled" is given, in
// which case eps_e = e * h_e * max(|u|+c) / P  (an O(h/P) viscosity).
// =============================================================================
#include "02_Basis1D.h"
#include "04_Mesh1D.h"
#include "05_CGAssembly.h"
#include "06_TimeIntegration.h"
#include "07_Output.h"
#include "09_CompressibleNS.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
double smoothStep(double x, double x0, double delta)
{
    if (delta <= 0.0)
        return x < x0 ? 0.0 : 1.0;
    return 0.5 * (1.0 + std::tanh((x - x0) / delta));
}
}   // namespace

int main(int argc, char **argv)
{
    Args args(argc, argv);
    if (args.has("--help"))
    {
        std::cout << "options: --case entropywave|sod|shuosher --nel N --P P --Q Q --T T --cfl c\n"
                     "         --mu mu --Pr Pr --gamma g --ic-smooth delta\n"
                     "         --svv --eps e [--scaled] --Mcut M --kernel exp|step --form book|kirby\n"
                     "         --svv-apply split|galerkin --nosvv-rho (do not filter density) --out tag --nprint n\n";
        return 0;
    }

    const std::string cas = args.getString("--case", "sod");
    const int nel         = args.getInt("--nel", 40);
    const int P           = args.getInt("--P", 8);
    const int Q           = args.getInt("--Q", P + 2);
    const double cfl      = args.getDouble("--cfl", 0.4);
    const int nprint      = args.getInt("--nprint", 200);
    const double delta    = args.getDouble("--ic-smooth", 0.0);

    GasProperties gas;
    gas.gamma = args.getDouble("--gamma", 1.4);
    gas.mu    = args.getDouble("--mu", 0.0);
    gas.Pr    = args.getDouble("--Pr", 0.72);

    CNSSVVOptions svv;
    svv.enabled   = args.has("--svv");
    svv.eps       = args.getDouble("--eps", 0.05);
    svv.scaled    = args.has("--scaled");
    svv.Mcut      = args.getInt("--Mcut", -1);
    svv.kernel    = parseKernelType(args.getString("--kernel", "exp"));
    svv.form      = parseSVVForm(args.getString("--form", "book"));
    svv.onDensity = !args.has("--nosvv-rho");
    svv.apply     = parseSVVApplication(args.getString("--svv-apply", "split"));

    // ---- case definition --------------------------------------------------------
    double xa, xb, T;
    BCType bc;
    std::function<double(double)> rho0, u0, p0;
    std::function<double(double, double)> rhoExact;   // (x,t) or empty

    if (cas == "entropywave")
    {
        xa = -1.0, xb = 1.0, T = 2.0, bc = BCType::Periodic;
        rho0     = [](double x) { return 1.0 + 0.2 * std::sin(M_PI * x); };
        u0       = [](double) { return 1.0; };
        p0       = [](double) { return 1.0; };
        rhoExact = [](double x, double t) { return 1.0 + 0.2 * std::sin(M_PI * (x - t)); };
    }
    else if (cas == "sod")
    {
        xa = 0.0, xb = 1.0, T = 0.2, bc = BCType::Dirichlet;
        rho0 = [delta](double x) { return 1.0 + (0.125 - 1.0) * smoothStep(x, 0.5, delta); };
        u0   = [](double) { return 0.0; };
        p0   = [delta](double x) { return 1.0 + (0.1 - 1.0) * smoothStep(x, 0.5, delta); };
    }
    else if (cas == "shuosher")
    {
        xa = -5.0, xb = 5.0, T = 1.8, bc = BCType::Dirichlet;
        rho0 = [delta](double x) {
            const double s = smoothStep(x, -4.0, delta);
            return 3.857143 * (1.0 - s) + (1.0 + 0.2 * std::sin(5.0 * x)) * s;
        };
        u0 = [delta](double x) { return 2.629369 * (1.0 - smoothStep(x, -4.0, delta)); };
        p0 = [delta](double x) { return 10.33333 * (1.0 - smoothStep(x, -4.0, delta)) + smoothStep(x, -4.0, delta); };
    }
    else
    {
        std::cerr << "unknown case " << cas << "\n";
        return 1;
    }
    T = args.getDouble("--T", T);

    std::string tag = args.getString("--out", cas + (svv.enabled ? "_svv" : "_nosvv"));

    // ---- discretisation ---------------------------------------------------------
    Mesh1D mesh = Mesh1D::uniform(xa, xb, nel);
    Basis1D basis(P, Q);
    CGAssembly A(mesh, basis, bc);
    CompressibleNS1D solver(A, gas, svv);

    Vector U = solver.projectPrimitive(rho0, u0, p0);
    if (bc == BCType::Dirichlet)
    {
        // pin the boundary dofs to the exact primitive states
        for (double x : {xa, xb})
        {
            const int d = (x == xa) ? A.dof(0, 0) : A.dof(nel - 1, P);
            const double r = rho0(x), uu = u0(x), pp = p0(x);
            solver.block(U, 0)[d] = r;
            solver.block(U, 1)[d] = r * uu;
            solver.block(U, 2)[d] = pp / (gas.gamma - 1.0) + 0.5 * r * uu * uu;
        }
    }

    FieldStats st0 = solver.stats(U);
    std::printf("Compressible %s 1D + SVV, case=%s\n", gas.mu > 0 ? "Navier-Stokes" : "Euler", cas.c_str());
    std::printf("  domain=[%g,%g]  elements=%d  P=%d  Q=%d  dofs/var=%d  T=%g\n", xa, xb, nel, P, Q, A.nGlobal, T);
    std::printf("  gas: gamma=%g  mu=%g  Pr=%g  kappa=%g\n", gas.gamma, gas.mu, gas.Pr, gas.kappa());
    std::printf("  SVV: %s", svv.enabled ? "on" : "off");
    if (svv.enabled)
    {
        std::printf("  eps=%g%s  M_SVV=%d  form=%s  apply=%s  density=%s  (eps_eff at t=0: %.4g)\n  kernel F_k:", svv.eps,
                    svv.scaled ? " (scaled: eps_e = eps*h*(|u|+c)/P)" : "", svv.Mcut < 0 ? P / 2 : svv.Mcut,
                    svv.form == SVVForm::Book ? "book" : "kirby", svv.apply == SVVApplication::Split ? "split" : "galerkin",
                    svv.onDensity ? "yes" : "no", solver.effectiveEps(st0));
        for (double f : solver.kernel())
            std::printf(" %.3f", f);
    }
    std::printf("\n  initial: min rho=%.4g  min p=%.4g  max(|u|+c)=%.4g\n", st0.minRho, st0.minP, st0.maxWaveSpeed);

    // ---- time loop ---------------------------------------------------------------
    RHSFunction R = [&](const Vector &u, double t) { return solver.rhs(u, t); };
    double t      = 0.0;
    int step      = 0;
    double dt     = solver.stableTimeStep(U, cfl, 0.8);
    std::printf("  initial dt=%.3e\n", dt);
    bool ok = true;
    try
    {
        while (t < T - 1e-14)
        {
            if (step % 10 == 0)
                dt = solver.stableTimeStep(U, cfl, 0.8);
            const double dtn = std::min(dt, T - t);
            rk4Step(U, t, dtn, R);
            solver.applySplitSVV(U, dtn);
            t += dtn;
            ++step;
            if (step % nprint == 0 || t >= T - 1e-14)
            {
                const FieldStats st = solver.stats(U);
                std::printf("  step %6d  t=%.4f  dt=%.2e  min rho=%.4f  min p=%.4f  max(|u|+c)=%.3f\n", step, t,
                            dtn, st.minRho, st.minP, st.maxWaveSpeed);
            }
        }
    }
    catch (const std::runtime_error &err)
    {
        std::printf("  ABORT at step %d: %s\n", step, err.what());
        ok = false;
    }

    // ---- output -------------------------------------------------------------------
    Vector x, rho, u, p;
    solver.samplePrimitive(U, 40, x, rho, u, p);
    const std::string csv = "results/" + tag + ".csv";
    writeCSV(csv, {{"x", &x}, {"rho", &rho}, {"u", &u}, {"p", &p}});
    Vector rhoHat(solver.block(U, 0), solver.block(U, 0) + A.nGlobal);
    writeSpectrumCSV("results/" + tag + "_spectrum.csv", A, rhoHat);

    std::printf("  final t=%.4f after %d steps; TV(rho)=%.4f  TV(u)=%.4f  TV(p)=%.4f\n", t, step, totalVariation(rho),
                totalVariation(u), totalVariation(p));
    if (rhoExact)
    {
        std::function<double(double)> ex = [&](double xx) { return rhoExact(xx, t); };
        std::printf("  L2 error of rho vs exact: %.6e\n", A.l2Norm(rhoHat, &ex));
    }
    std::printf("  wrote %s and results/%s_spectrum.csv\n", csv.c_str(), tag.c_str());
    return ok ? 0 : 2;
}
