#include "09_CompressibleNS.h"
#include "07_Output.h"

#include <algorithm>

#include <cmath>
#include <sstream>
#include <stdexcept>

namespace
{
Matrix buildMassBC(const CGAssembly &A)
{
    Matrix M = A.massMatrix();
    A.applyDirichletRows(M);
    return M;
}
}   // namespace

CompressibleNS1D::CompressibleNS1D(const CGAssembly &A, const GasProperties &gas, const CNSSVVOptions &svv)
    : A_(A), gas_(gas), svv_(svv), kernel_(svvKernel(A.P, svv.Mcut < 0 ? A.P / 2 : svv.Mcut, svv.kernel)),
      Lref_(svvOperator(A.basis, kernel_, svv.form)), massLU_(buildMassBC(A)), bubbleFilter_(A.basis, Lref_)
{
    if (svv_.Mcut < 0)
        svv_.Mcut = A.P / 2;

    // Spectral radii for the time-step estimate (linear operators, computed once)
    if (svv_.enabled && svv_.apply == SVVApplication::Galerkin)
        rhoSVV_ = A.spectralRadius(A.operatorMatrix(Lref_));
    if (gas_.mu > 0.0)
        rhoLap_ = A.spectralRadius(A.laplacianMatrix());
}

Vector CompressibleNS1D::projectPrimitive(const std::function<double(double)> &rho,
                                          const std::function<double(double)> &u,
                                          const std::function<double(double)> &p) const
{
    const double g = gas_.gamma;
    auto m = [&](double x) { return rho(x) * u(x); };
    auto E = [&](double x) { return p(x) / (g - 1.0) + 0.5 * rho(x) * u(x) * u(x); };

    Vector U(size());
    const Vector r = A_.l2Projection(rho), mm = A_.l2Projection(m), ee = A_.l2Projection(E);
    std::copy(r.begin(), r.end(), block(U, 0));
    std::copy(mm.begin(), mm.end(), block(U, 1));
    std::copy(ee.begin(), ee.end(), block(U, 2));
    return U;
}

CompressibleNS1D::PointState CompressibleNS1D::pointState(double rho, double m, double E, double drho,
                                                          double dm, double dE, double x, double t) const
{
    const double g = gas_.gamma;
    PointState s;
    s.rho = rho;
    s.u   = m / rho;
    s.p   = (g - 1.0) * (E - 0.5 * m * s.u);
    if (!(rho > 0.0) || !(s.p > 0.0))
    {
        std::ostringstream os;
        os << "negative density/pressure at x=" << x << ", t=" << t << " (rho=" << rho << ", p=" << s.p
           << "). The solution blew up: add SVV / viscosity, lower P, or smooth the IC.";
        throw std::runtime_error(os.str());
    }
    s.T = s.p / (rho * gas_.R);
    // chain rule from the exact polynomial gradients of the conserved variables
    s.ux            = (dm - s.u * drho) / rho;
    const double px = (g - 1.0) * (dE - s.u * dm + 0.5 * s.u * s.u * drho);
    s.Tx            = (px - gas_.R * s.T * drho) / (rho * gas_.R);
    return s;
}

Vector CompressibleNS1D::rhs(const Vector &U, double t) const
{
    const Basis1D &b = A_.basis;
    const int Q      = b.Q;
    const int n      = A_.nGlobal;

    Vector r(size(), 0.0);
    Vector F1(Q), F2(Q), F3(Q), G2(Q), G3(Q);

    // SVV amplitude (may depend on the state when "scaled")
    double epsOverJFactor = 0.0;   // multiplies L_ref; per element divided by J or not
    FieldStats st;
    if (svv_.enabled)
    {
        if (svv_.scaled)
            st = stats(U);
    }

    for (int e = 0; e < A_.nel; ++e)
    {
        const double J    = A_.mesh.J(e);
        const Vector rhoh = A_.gather(block(U, 0), e);
        const Vector mh   = A_.gather(block(U, 1), e);
        const Vector Eh   = A_.gather(block(U, 2), e);

        const Vector rho = b.evaluate(rhoh), m = b.evaluate(mh), E = b.evaluate(Eh);
        // d/dx = (1/J) d/dxi
        const Vector drho = scaled(b.evaluateDeriv(rhoh), 1.0 / J);
        const Vector dm   = scaled(b.evaluateDeriv(mh), 1.0 / J);
        const Vector dE   = scaled(b.evaluateDeriv(Eh), 1.0 / J);

        for (int q = 0; q < Q; ++q)
        {
            const PointState s =
                pointState(rho[q], m[q], E[q], drho[q], dm[q], dE[q], A_.mesh.x(e, b.zq[q]), t);
            F1[q] = m[q];
            F2[q] = m[q] * s.u + s.p;
            F3[q] = s.u * (E[q] + s.p);

            const double tau = (4.0 / 3.0) * gas_.mu * s.ux;
            const double qh  = -gas_.kappa() * s.Tx;
            G2[q]            = tau;
            G3[q]            = s.u * tau - qh;
        }

        // convection: -(psi_i, dF/dx) = -S f^
        Vector r1 = scaled(b.S * b.project(F1), -1.0);
        Vector r2 = scaled(b.S * b.project(F2), -1.0);
        Vector r3 = scaled(b.S * b.project(F3), -1.0);

        // viscous: -(psi_i', G)
        if (gas_.mu > 0.0)
        {
            axpy(-1.0, b.innerProductDeriv(G2), r2);
            axpy(-1.0, b.innerProductDeriv(G3), r3);
        }

        // SVV (Galerkin application): -(eps_e / J) L_ref u^
        if (svv_.enabled && svv_.apply == SVVApplication::Galerkin)
        {
            if (svv_.scaled)
                epsOverJFactor = 2.0 * svv_.eps * st.maxWaveSpeed / A_.P;   // eps_e = mu h c/P
            else
                epsOverJFactor = svv_.eps / J;
            if (svv_.onDensity)
                axpy(-epsOverJFactor, Lref_ * rhoh, r1);
            axpy(-epsOverJFactor, Lref_ * mh, r2);
            axpy(-epsOverJFactor, Lref_ * Eh, r3);
        }

        A_.scatterAdd(r1, e, r.data());
        A_.scatterAdd(r2, e, r.data() + n);
        A_.scatterAdd(r3, e, r.data() + 2 * n);
    }

    Vector dU(size());
    for (int k = 0; k < nvar; ++k)
    {
        Vector rk(r.begin() + k * n, r.begin() + (k + 1) * n);
        A_.zeroDirichlet(rk);
        const Vector dk = massLU_.solve(rk);
        std::copy(dk.begin(), dk.end(), block(dU, k));
    }
    return dU;
}

void CompressibleNS1D::applySplitSVV(Vector &U, double dt)
{
    if (!svv_.enabled || svv_.apply != SVVApplication::Split)
        return;
    const double cmax = svv_.scaled ? stats(U).maxWaveSpeed : 0.0;
    for (int e = 0; e < A_.nel; ++e)
    {
        // element system (J M) dU/dt = -(eps_e/J) L U, divided by J -> c = dt eps_e / J^2
        // scaled: eps_e / J = 2 mu_svv c_max / P, so c = dt * 2 mu_svv c_max / (P J)
        const double J = A_.mesh.J(e);
        const double c = svv_.scaled ? dt * 2.0 * svv_.eps * cmax / (A_.P * J) : dt * svv_.eps / (J * J);
        for (int k = 0; k < nvar; ++k)
        {
            if (k == 0 && !svv_.onDensity)
                continue;
            double *g       = block(U, k);
            const Vector ue = bubbleFilter_.apply(A_.gather(g, e), c);
            for (int p = 1; p < A_.P; ++p)
                g[A_.dof(e, p)] = ue[p];
        }
    }
}

FieldStats CompressibleNS1D::stats(const Vector &U) const
{
    const Basis1D &b = A_.basis;
    FieldStats st;
    st.minRho = st.minP = 1e300;
    for (int e = 0; e < A_.nel; ++e)
    {
        const Vector rho = b.evaluate(A_.gather(block(U, 0), e));
        const Vector m   = b.evaluate(A_.gather(block(U, 1), e));
        const Vector E   = b.evaluate(A_.gather(block(U, 2), e));
        for (int q = 0; q < b.Q; ++q)
        {
            const double u = m[q] / rho[q];
            const double p = (gas_.gamma - 1.0) * (E[q] - 0.5 * m[q] * u);
            st.minRho      = std::min(st.minRho, rho[q]);
            st.minP        = std::min(st.minP, p);
            if (rho[q] > 0.0 && p > 0.0)
            {
                const double c   = std::sqrt(gas_.gamma * p / rho[q]);
                st.maxWaveSpeed  = std::max(st.maxWaveSpeed, std::fabs(u) + c);
                const double nu  = (4.0 / 3.0) * gas_.mu / rho[q];
                const double alp = gas_.kappa() / (rho[q] * gas_.cv());
                st.maxDiffusivity = std::max(st.maxDiffusivity, std::max(nu, alp));
            }
        }
    }
    return st;
}

double CompressibleNS1D::effectiveEps(const FieldStats &st) const
{
    if (!svv_.enabled)
        return 0.0;
    if (!svv_.scaled)
        return svv_.eps;
    // representative value for reporting: use the first element's h
    return svv_.eps * A_.mesh.h(0) * st.maxWaveSpeed / A_.P;
}

double CompressibleNS1D::stableTimeStep(const Vector &U, double cflAdv, double safetyDiff) const
{
    const FieldStats st = stats(U);
    const double hmin   = A_.minSpacing();
    double dt           = cflAdv * hmin / std::max(st.maxWaveSpeed, 1e-300);

    // diffusive limit: RK4 real-axis stability ~ 2.785; use "safetyDiff * 2.785"
    double lambda = 0.0;
    if (gas_.mu > 0.0)
        lambda += st.maxDiffusivity * rhoLap_;
    if (svv_.enabled && svv_.apply == SVVApplication::Galerkin)
    {
        // effective eps/J is uniform for scaled mode; for const mode rhoSVV_ already
        // contains the 1/J factors, so multiply by eps.
        if (svv_.scaled)
            lambda += (2.0 * svv_.eps * st.maxWaveSpeed / A_.P) * rhoSVV_ * A_.mesh.J(0);
        else
            lambda += svv_.eps * rhoSVV_;
    }
    if (lambda > 0.0)
        dt = std::min(dt, safetyDiff * 2.785 / lambda);
    return dt;
}

void CompressibleNS1D::samplePrimitive(const Vector &U, int nsPerElem, Vector &x, Vector &rho, Vector &u,
                                       Vector &p) const
{
    Samples sr = sampleField(A_, block(U, 0), nsPerElem);
    Samples sm = sampleField(A_, block(U, 1), nsPerElem);
    Samples sE = sampleField(A_, block(U, 2), nsPerElem);
    x          = sr.x;
    rho        = sr.u;
    u.resize(x.size());
    p.resize(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        u[i] = sm.u[i] / sr.u[i];
        p[i] = (gas_.gamma - 1.0) * (sE.u[i] - 0.5 * sm.u[i] * u[i]);
    }
}
