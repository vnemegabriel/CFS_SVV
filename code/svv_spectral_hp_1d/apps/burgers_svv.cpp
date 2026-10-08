// =============================================================================
// burgers_svv  --  reproduces Karniadakis & Sherwin Fig. 6.27 (p. 351):
//
//   "Solution of the inviscid Burgers equation at time T = 0.5 using continuous
//    Galerkin (a) without, and (b) with SVV. Five equally-spaced elements
//    spanning [-1,1] were used, each of which contained sixteen modes.
//    A wave cut-off of M_SVV = 8 and amplitude eps = 1/16 were applied."
//
// Usage examples
//   ./burgers_svv                      # with SVV (book parameters)
//   ./burgers_svv --nosvv              # without SVV -> Fig. 6.27(a)
//   ./burgers_svv --Mcut 12 --eps 0.02 # play with the kernel
//   ./burgers_svv --form kirby         # eq. (6.5.17) variant
//   ./burgers_svv --Q 24               # over-integrate the u^2 term
//
// The initial condition is u0 = -sin(pi x) (shock forms at x = 0 at t = 1/pi,
// so at T = 0.5 there is a standing shock at the origin).
// =============================================================================
#include "02_Basis1D.h"
#include "04_Mesh1D.h"
#include "05_CGAssembly.h"
#include "06_TimeIntegration.h"
#include "07_Output.h"
#include "08_Burgers.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char **argv)
{
    Args args(argc, argv);
    if (args.has("--help"))
    {
        std::cout << "options: --nel N --P P --Q Q --T T --cfl c --nosvv --eps e --Mcut M\n"
                     "         --kernel exp|step --form book|kirby --conv cons|nc|skew\n"
                     "         --svv-apply split|galerkin --ic sine|step --dt dt --out tag\n";
        return 0;
    }

    // ---- parameters (defaults = the book's test) ------------------------------
    const int nel      = args.getInt("--nel", 5);
    const int P        = args.getInt("--P", 15);          // 16 modes
    const int Q        = args.getInt("--Q", P + 1);       // GLL points; P+1 = classic collocation choice
    const double T     = args.getDouble("--T", 0.5);
    const double cfl   = args.getDouble("--cfl", 0.5);
    const int nprint   = args.getInt("--nprint", 500);
    const std::string ic = args.getString("--ic", "sine");

    BurgersOptions opt;
    opt.useSVV = !args.has("--nosvv");
    opt.eps    = args.getDouble("--eps", 1.0 / 16.0);
    opt.Mcut   = args.getInt("--Mcut", 8);
    opt.kernel = parseKernelType(args.getString("--kernel", "exp"));
    opt.form   = parseSVVForm(args.getString("--form", "book"));
    {
        const std::string c = args.getString("--conv", "cons");
        opt.conv = c == "skew" ? ConvectiveForm::Skew : c == "nc" ? ConvectiveForm::NonConservative : ConvectiveForm::Conservative;
        opt.apply = parseSVVApplication(args.getString("--svv-apply", "split"));
    }

    const std::string tag = args.getString(
        "--out", !opt.useSVV ? "burgers_nosvv" : opt.apply == SVVApplication::Split ? "burgers_svv" : "burgers_svv_galerkin");

    // ---- discretisation -------------------------------------------------------
    Mesh1D mesh = Mesh1D::uniform(-1.0, 1.0, nel);
    Basis1D basis(P, Q);
    CGAssembly A(mesh, basis, BCType::Dirichlet);
    BurgersSolver solver(A, opt);

    auto u0 = [&](double x) {
        if (ic == "step")
            return x < 0.0 ? 1.0 : 0.0;
        return -std::sin(M_PI * x);
    };
    Vector uhat = A.l2Projection(u0);
    uhat[A.dof(0, 0)]       = u0(mesh.xmin());   // enforce Dirichlet values exactly
    uhat[A.dof(nel - 1, P)] = u0(mesh.xmax());

    // ---- time step ------------------------------------------------------------
    const double hmin   = A.minSpacing();
    const double umax   = std::max(maxAbs(sampleField(A, uhat).u), 1e-12);
    const double dtAdv  = cfl * hmin / (1.2 * umax);
    const double rhoSVV = opt.apply == SVVApplication::Galerkin ? solver.svvSpectralRadius() : 0.0;
    const double dtSVV  = rhoSVV > 0.0 ? 0.8 * 2.785 / rhoSVV : 1e300;
    double dt           = std::min(dtAdv, dtSVV);
    if (args.has("--dt"))
        dt = args.getDouble("--dt", dt);   // manual override (stability experiments)
    const int nsteps    = static_cast<int>(std::ceil(T / dt));
    dt                  = T / nsteps;

    std::printf("Burgers + SVV (continuous Galerkin spectral/hp)\n");
    std::printf("  elements=%d  P=%d (modes=%d)  Q=%d  global dofs=%d\n", nel, P, P + 1, Q, A.nGlobal);
    std::printf("  convective form: %s\n", opt.conv == ConvectiveForm::Skew ? "skew-symmetric"
                                            : opt.conv == ConvectiveForm::NonConservative ? "non-conservative"
                                                                                          : "conservative");
    std::printf("  SVV: %s  eps=%.5g  M_SVV=%d  kernel=%s  form=%s  application=%s\n", opt.useSVV ? "on" : "off",
                opt.eps, opt.Mcut, opt.kernel == SVVKernelType::Exponential ? "exp" : "step",
                opt.form == SVVForm::Book ? "book" : "kirby",
                opt.apply == SVVApplication::Split ? "split (element-wise implicit bubble filter)"
                                                   : "Galerkin (weak-form residual, global mass matrix)");
    if (opt.useSVV)
    {
        std::printf("  kernel F_k:");
        for (double f : solver.kernel())
            std::printf(" %.3f", f);
        std::printf("\n");
    }
    std::printf("  h_min(GLL)=%.4g  dt_adv=%.3e  dt_svv=%.3e (rho=%.3e)  -> dt=%.3e, %d steps\n", hmin,
                dtAdv, dtSVV, rhoSVV, dt, nsteps);

    // ---- time loop ------------------------------------------------------------
    RHSFunction R = [&](const Vector &u, double t) { return solver.rhs(u, t); };
    double t      = 0.0;
    for (int n = 0; n < nsteps; ++n)
    {
        rk4Step(uhat, t, dt, R);
        solver.applySplitSVV(uhat, dt);   // no-op unless --svv-apply split
        t += dt;
        if (!std::isfinite(maxAbs(uhat)) || maxAbs(uhat) > 1e3)
        {
            std::printf("  BLOW-UP at t=%.4f (step %d)\n", t, n + 1);
            break;
        }
        if ((n + 1) % nprint == 0 || n + 1 == nsteps)
            std::printf("  step %6d  t=%.4f  ||u||_L2=%.6f  max|u|=%.4f\n", n + 1, t, A.l2Norm(uhat),
                        maxAbs(sampleField(A, uhat).u));
    }

    // ---- diagnostics + output ---------------------------------------------------
    Samples s = sampleField(A, uhat, 80);
    // exact TV: step IC -> 1; sine IC -> 4 for all t (0 -> 1 -> -1 -> 0; for t > 1/pi the
    // characteristic from x0 = -1/2 reaches the shock at t = 1/2, so u(0-) = 1 at T = 0.5)
    std::printf("  total variation of u_h at T: %.4f  (exact solution: %.4f)\n", totalVariation(s.u),
                ic == "step" ? 1.0 : 4.0);

    const std::string csv = "results/" + tag + ".csv";
    writeCSV(csv, {{"x", &s.x}, {"u", &s.u}});
    writeSpectrumCSV("results/" + tag + "_spectrum.csv", A, uhat);
    std::printf("  wrote %s and results/%s_spectrum.csv\n", csv.c_str(), tag.c_str());
    return 0;
}
