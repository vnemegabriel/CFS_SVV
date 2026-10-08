// =============================================================================
// 09_CompressibleNS.h  --  1D compressible Navier-Stokes with SVV
//
// Conserved variables U = (rho, m = rho u, E), ideal gas p = (g-1)(E - m^2/(2 rho)).
//
//      U_t + F(U)_x = G(U, U_x)_x + eps d/dx [ F * U_x ]
//
//      F = ( m,  m u + p,  u (E + p) )                        inviscid flux
//      G = ( 0,  tau,      u tau - q ),  tau = 4/3 mu u_x,  q = -kappa T_x
//
// with kappa = mu c_p / Pr and T = p / (rho R).  Setting mu = 0 gives Euler.
//
// Weak form (test function v, each conserved variable separately):
//
//      (v, U_t) = -(v, F_x) - (v_x, G) - eps (v_x, F * U_x)
//
// so the C0 Galerkin RHS per element is
//      -S f^        (f^ = projection of F(U_q))
//      -dB^T W G_q  (viscous, one d/dx and one Jacobian cancel)
//      -(eps/J) L_svv u^
//
// GRADIENTS. U_x is computed exactly from the modal expansion (U is a polynomial
// in each element); u_x, p_x, T_x follow by the chain rule at the quadrature
// points. This avoids differentiating non-polynomial quantities.
//
// SVV AMPLITUDE. Two options:
//   * constant eps, as in the Burgers test;
//   * "scaled": eps_e = mu_svv * h_e * max(|u|+c) / P   -- i.e. an O(h/P)
//     viscosity, the natural scaling for spectral-element under-resolution
//     (Kirby & Karniadakis 2002, Nektar++'s SVVDiffCoeff). Because the SVV
//     matrix carries a 1/J_e = 2/h_e factor, eps_e/J_e = 2 mu_svv c_max / P is
//     independent of the element size. The damping RATE is not: dividing by
//     the element mass J_e M gives eps_e/J_e^2 ~ 1/h (O(h/P) viscosity acting
//     on modes of wavenumber ~P/h).
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"
#include "03_SVV.h"
#include "05_CGAssembly.h"

#include <functional>

struct GasProperties
{
    double gamma = 1.4;
    double R     = 1.0;     // gas constant (non-dimensional)
    double mu    = 0.0;     // dynamic viscosity (0 -> Euler)
    double Pr    = 0.72;
    double cv() const { return R / (gamma - 1.0); }
    double cp() const { return gamma * R / (gamma - 1.0); }
    double kappa() const { return mu * cp() / Pr; }
};

struct CNSSVVOptions
{
    bool enabled   = false;
    double eps     = 0.0;      // constant eps, or mu_svv if scaled
    bool scaled    = false;
    int Mcut       = -1;       // -1 -> P/2
    SVVKernelType kernel = SVVKernelType::Exponential;
    SVVForm form   = SVVForm::Book;
    SVVApplication apply = SVVApplication::Split;
    bool onDensity = true;     // apply SVV to rho as well as to m and E
};

struct FieldStats
{
    double maxWaveSpeed = 0.0;   // max(|u| + c)
    double minRho       = 0.0;
    double minP         = 0.0;
    double maxDiffusivity = 0.0; // max over points of max(4/3 mu/rho, kappa/(rho cv))
};

class CompressibleNS1D
{
public:
    CompressibleNS1D(const CGAssembly &A, const GasProperties &gas, const CNSSVVOptions &svv);

    static constexpr int nvar = 3;
    int nGlobal() const { return A_.nGlobal; }
    int size() const { return nvar * A_.nGlobal; }

    // Block layout of the global state: U = [rho^ | m^ | E^]
    const double *block(const Vector &U, int k) const { return U.data() + k * A_.nGlobal; }
    double *block(Vector &U, int k) const { return U.data() + k * A_.nGlobal; }

    // Project primitive initial condition (rho, u, p)(x) onto the C0 space.
    Vector projectPrimitive(const std::function<double(double)> &rho,
                            const std::function<double(double)> &u,
                            const std::function<double(double)> &p) const;

    // dU/dt = M_bc^{-1} r(U)
    Vector rhs(const Vector &U, double t) const;

    // Split application of SVV: call once after every time step of size dt.
    void applySplitSVV(Vector &U, double dt);

    FieldStats stats(const Vector &U) const;

    // dt from advective CFL and from the (physical + SVV) diffusive limits
    double stableTimeStep(const Vector &U, double cflAdv, double safetyDiff) const;

    // effective SVV coefficient used at the given state (const or scaled)
    double effectiveEps(const FieldStats &st) const;

    const Vector &kernel() const { return kernel_; }
    const GasProperties &gas() const { return gas_; }

    // Evaluate primitive variables of the discrete solution at sample points.
    void samplePrimitive(const Vector &U, int nsPerElem, Vector &x, Vector &rho, Vector &u, Vector &p) const;

private:
    struct PointState
    {
        double rho, u, p, T, ux, Tx;
    };
    PointState pointState(double rho, double m, double E, double drho, double dm, double dE,
                          double x, double t) const;

    const CGAssembly &A_;
    GasProperties gas_;
    CNSSVVOptions svv_;
    Vector kernel_;
    Matrix Lref_;
    LUSolver massLU_;
    SVVBubbleFilter bubbleFilter_;
    double rhoSVV_ = 0.0;   // spectral radius of M^{-1} L_svv (unit eps)
    double rhoLap_ = 0.0;   // spectral radius of M^{-1} Laplacian
};
