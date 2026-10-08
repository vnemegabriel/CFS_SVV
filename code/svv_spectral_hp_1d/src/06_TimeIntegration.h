// =============================================================================
// 06_TimeIntegration.h  --  explicit Runge-Kutta for  du/dt = R(u, t)
//
// The spatial discretisation produces  M du^/dt = r(u^)  and the "RHS function"
// handed to the integrator already includes the mass-matrix solve, i.e.
// R(u^) = M^{-1} r(u^).  This "method of lines" split is the same one Nektar++
// uses (TimeIntegrationScheme acting on DoOdeRhs + DoOdeProjection).
//
// Stability reminder for the linear test equation y' = lambda y:
//   * classical RK4 is stable for real lambda*dt in (-2.785, 0]
//   * SSP-RK3 (Shu-Osher) is stable for real lambda*dt in (-2.51, 0]
// The SVV term is a (filtered) Laplacian, so its most negative eigenvalue
// sets a diffusive dt limit  dt < 2.785 / (eps * rho(M^{-1} L_svv)).
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"

#include <functional>

using RHSFunction = std::function<Vector(const Vector &, double)>;

inline void rk4Step(Vector &u, double t, double dt, const RHSFunction &R)
{
    const Vector k1 = R(u, t);
    Vector tmp(u);
    axpy(0.5 * dt, k1, tmp);
    const Vector k2 = R(tmp, t + 0.5 * dt);
    tmp             = u;
    axpy(0.5 * dt, k2, tmp);
    const Vector k3 = R(tmp, t + 0.5 * dt);
    tmp             = u;
    axpy(dt, k3, tmp);
    const Vector k4 = R(tmp, t + dt);
    for (size_t i = 0; i < u.size(); ++i)
        u[i] += dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

inline void ssprk3Step(Vector &u, double t, double dt, const RHSFunction &R)
{
    Vector u1(u);
    axpy(dt, R(u, t), u1);
    Vector u2(u1);
    axpy(dt, R(u1, t + dt), u2);
    for (size_t i = 0; i < u.size(); ++i)
        u2[i] = 0.75 * u[i] + 0.25 * u2[i];
    Vector u3(u2);
    axpy(dt, R(u2, t + 0.5 * dt), u3);
    for (size_t i = 0; i < u.size(); ++i)
        u[i] = (1.0 / 3.0) * u[i] + (2.0 / 3.0) * u3[i];
}
