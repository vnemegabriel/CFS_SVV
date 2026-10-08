// cfs1d.cpp — 1D Euler, weak DG, same call sequence as Nektar++ CompressibleFlowSolver
// (EulerCFE, AdvectionWeakDG, LaxFriedrichs or Roe Riemann solver, explicit RK4).
// Names follow Nektar++ functions so the trace can be compared line by line.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <string>
#include <vector>

using Vec = std::vector<double>;
using Mat = std::vector<Vec>;
const double PI = std::acos(-1.0), GAMMA = 1.4;

// ---------------- LibUtilities::Foundations (Points, Basis) ----------------
double Jacobi(int n, double a, double b, double x)            // P_n^{a,b}(x)
{
    if (n == 0) return 1.0;
    double p0 = 1.0, p1 = 0.5 * (a - b + (a + b + 2) * x);
    for (int k = 2; k <= n; ++k)
    {
        double c = 2 * k + a + b;
        double a1 = 2 * k * (k + a + b) * (c - 2);
        double a2 = (c - 1) * (a * a - b * b);
        double a3 = (c - 2) * (c - 1) * c;
        double a4 = 2 * (k + a - 1) * (k + b - 1) * c;
        double p2 = ((a2 + a3 * x) * p1 - a4 * p0) / a1;
        p0 = p1; p1 = p2;
    }
    return p1;
}
double dJacobi(int n, double a, double b, double x)
{
    return n == 0 ? 0.0 : 0.5 * (n + a + b + 1) * Jacobi(n - 1, a + 1, b + 1, x);
}

// Gauss-Lobatto-Legendre points and weights (eGaussLobattoLegendre), Q points
void GLL(int Q, Vec &z, Vec &w)
{
    int N = Q - 1;
    z.assign(Q, 0); w.assign(Q, 0);
    for (int i = 0; i < Q; ++i)
    {
        double x = -std::cos(PI * i / N);
        for (int it = 0; it < 100; ++it)
        {
            double PN = Jacobi(N, 0, 0, x), PNm1 = Jacobi(N - 1, 0, 0, x);
            double dx = (x * PN - PNm1) / (Q * PN);
            x -= dx;
            if (std::abs(dx) < 1e-16) break;
        }
        z[i] = x;
        double PN = Jacobi(N, 0, 0, x);
        w[i] = 2.0 / (N * (N + 1) * PN * PN);
    }
}

// Modified_A basis: phi_0 = (1-x)/2, phi_1 = (1+x)/2, phi_p = (1-x)(1+x)/4 P^{1,1}_{p-2}
double phi(int p, double x)
{
    if (p == 0) return 0.5 * (1 - x);
    if (p == 1) return 0.5 * (1 + x);
    return 0.25 * (1 - x) * (1 + x) * Jacobi(p - 2, 1, 1, x);
}
double dphi(int p, double x)
{
    if (p == 0) return -0.5;
    if (p == 1) return 0.5;
    return -0.5 * x * Jacobi(p - 2, 1, 1, x) + 0.25 * (1 - x * x) * dJacobi(p - 2, 1, 1, x);
}

Mat inverse(Mat M)
{
    int N = M.size();
    Mat I(N, Vec(N, 0));
    for (int i = 0; i < N; ++i) I[i][i] = 1;
    for (int k = 0; k < N; ++k)
    {
        int p = k;
        for (int r = k + 1; r < N; ++r) if (std::abs(M[r][k]) > std::abs(M[p][k])) p = r;
        std::swap(M[k], M[p]); std::swap(I[k], I[p]);
        double d = M[k][k];
        for (int j = 0; j < N; ++j) { M[k][j] /= d; I[k][j] /= d; }
        for (int r = 0; r < N; ++r)
        {
            if (r == k) continue;
            double m = M[r][k];
            for (int j = 0; j < N; ++j) { M[r][j] -= m * M[k][j]; I[r][j] -= m * I[k][j]; }
        }
    }
    return I;
}

// ---------------- StdRegions / LocalRegions: one segment ----------------
struct SegExp
{
    int ncoeffs, nq;
    double J;                 // dx/dxi = h/2 (SpatialDomains::GeomFactors)
    Vec z, w;
    Mat B, DB, Minv;          // B[q][p] = phi_p(z_q)

    SegExp(int P, int Q, double h) : ncoeffs(P + 1), nq(Q), J(0.5 * h)
    {
        GLL(Q, z, w);
        B.assign(Q, Vec(ncoeffs)); DB = B;
        for (int q = 0; q < Q; ++q)
            for (int p = 0; p < ncoeffs; ++p) { B[q][p] = phi(p, z[q]); DB[q][p] = dphi(p, z[q]); }
        Mat M(ncoeffs, Vec(ncoeffs, 0));                       // M = B^T Lambda(wJ) B
        for (int p = 0; p < ncoeffs; ++p)
            for (int r = 0; r < ncoeffs; ++r)
                for (int q = 0; q < Q; ++q) M[p][r] += B[q][p] * w[q] * J * B[q][r];
        Minv = inverse(M);
    }
    Vec BwdTrans(const Vec &c) const
    {
        Vec u(nq, 0);
        for (int q = 0; q < nq; ++q) for (int p = 0; p < ncoeffs; ++p) u[q] += B[q][p] * c[p];
        return u;
    }
    Vec IProductWRTBase(const Vec &f) const                   // B^T Lambda(wJ) f
    {
        Vec r(ncoeffs, 0);
        for (int p = 0; p < ncoeffs; ++p) for (int q = 0; q < nq; ++q) r[p] += B[q][p] * w[q] * J * f[q];
        return r;
    }
    Vec IProductWRTDerivBase(const Vec &f) const              // (dphi/dx, f) = DB^T Lambda(w) f
    {
        Vec r(ncoeffs, 0);
        for (int p = 0; p < ncoeffs; ++p) for (int q = 0; q < nq; ++q) r[p] += DB[q][p] * w[q] * f[q];
        return r;
    }
    Vec MultiplyByElmtInvMass(const Vec &c) const
    {
        Vec r(ncoeffs, 0);
        for (int p = 0; p < ncoeffs; ++p) for (int s = 0; s < ncoeffs; ++s) r[p] += Minv[p][s] * c[s];
        return r;
    }
    Vec FwdTrans(const Vec &u) const { return MultiplyByElmtInvMass(IProductWRTBase(u)); }
};

// ---------------- MultiRegions::DisContField (1D, uniform mesh) ----------------
struct DisContField
{
    int ne; bool periodic; double x0, h;
    SegExp exp;
    DisContField(int ne_, int P, double a, double b, bool per)
        : ne(ne_), periodic(per), x0(a), h((b - a) / ne_), exp(P, P + 2, (b - a) / ne_) {}
    int nTrace() const { return periodic ? ne : ne + 1; }
    double xq(int e, int q) const { return x0 + h * e + exp.J * (exp.z[q] + 1); }

    // trace k sits at x0 + k h; Fwd = left element (normal +x points out of it), Bwd = right element
    void GetFwdBwdTracePhys(const std::vector<Vec> &phys, Vec &Fwd, Vec &Bwd) const
    {
        int nt = nTrace(), Q = exp.nq;
        Fwd.assign(nt, 0); Bwd.assign(nt, 0);
        for (int k = 0; k < nt; ++k)
        {
            int eL = k - 1, eR = k;
            if (periodic) { eL = (k - 1 + ne) % ne; eR = k % ne; }
            if (eL < 0) eL = -1;
            if (eR >= ne) eR = -1;
            Fwd[k] = eL >= 0 ? phys[eL][Q - 1] : phys[eR][0];   // outflow (copy) boundary
            Bwd[k] = eR >= 0 ? phys[eR][0] : phys[eL][Q - 1];
        }
    }
    // adds the trace integral <phi, Fn n>: +phi(+1)Fn at right end, -phi(-1)Fn at left end
    void AddTraceIntegral(const Vec &Fn, std::vector<Vec> &coeffs) const
    {
        int nt = nTrace();
        for (int e = 0; e < ne; ++e)
        {
            int kL = e, kR = periodic ? (e + 1) % nt : e + 1;
            coeffs[e][1] += Fn[kR];   // Modified_A: only phi_1 is nonzero at xi = +1
            coeffs[e][0] -= Fn[kL];   // only phi_0 is nonzero at xi = -1
        }
    }
};

// ---------------- CompressibleFlowSolver (EulerCFE) ----------------
using Field = std::vector<std::vector<Vec>>;     // [variable][element][point or mode]
struct EulerCFE
{
    DisContField f;
    std::string upwind;
    double svvCoeff = 0, svvCutoff = 0.5;
    bool trace = false;

    EulerCFE(int ne, int P, double a, double b, bool per, std::string up)
        : f(ne, P, a, b, per), upwind(up) {}

    static double pressure(double r, double ru, double E) { return (GAMMA - 1) * (E - 0.5 * ru * ru / r); }

    // GetFluxVector: F(U) at every quadrature point
    Field GetFluxVector(const Field &U) const
    {
        Field F = U;
        for (int e = 0; e < f.ne; ++e)
            for (int q = 0; q < f.exp.nq; ++q)
            {
                double r = U[0][e][q], ru = U[1][e][q], E = U[2][e][q], u = ru / r, p = pressure(r, ru, E);
                F[0][e][q] = ru; F[1][e][q] = ru * u + p; F[2][e][q] = u * (E + p);
            }
        return F;
    }
    // RiemannSolver::Solve on the trace (1D: normal = +x, no rotation needed)
    void RiemannSolve(const std::vector<Vec> &Fwd, const std::vector<Vec> &Bwd, std::vector<Vec> &Fn) const
    {
        int nt = Fwd[0].size();
        Fn.assign(3, Vec(nt, 0));
        for (int k = 0; k < nt; ++k)
        {
            double L[3] = {Fwd[0][k], Fwd[1][k], Fwd[2][k]}, R[3] = {Bwd[0][k], Bwd[1][k], Bwd[2][k]};
            double uL = L[1] / L[0], pL = pressure(L[0], L[1], L[2]), cL = std::sqrt(GAMMA * pL / L[0]);
            double uR = R[1] / R[0], pR = pressure(R[0], R[1], R[2]), cR = std::sqrt(GAMMA * pR / R[0]);
            double FL[3] = {L[1], L[1] * uL + pL, uL * (L[2] + pL)};
            double FR[3] = {R[1], R[1] * uR + pR, uR * (R[2] + pR)};
            if (upwind == "LaxFriedrichs")
            {
                double lam = std::max(std::abs(uL) + cL, std::abs(uR) + cR);
                for (int i = 0; i < 3; ++i) Fn[i][k] = 0.5 * (FL[i] + FR[i]) - 0.5 * lam * (R[i] - L[i]);
                continue;
            }
            // Roe: F = (FL+FR)/2 - 1/2 sum |lambda_j| alpha_j r_j   (K&S App. E characteristic split)
            double sL = std::sqrt(L[0]), sR = std::sqrt(R[0]);
            double u = (sL * uL + sR * uR) / (sL + sR);
            double H = (sL * (L[2] + pL) / L[0] + sR * (R[2] + pR) / R[0]) / (sL + sR);
            double c = std::sqrt((GAMMA - 1) * (H - 0.5 * u * u));
            double dr = R[0] - L[0], du = uR - uL, dp = pR - pL, rho = sL * sR;
            double a1 = (dp - rho * c * du) / (2 * c * c), a2 = dr - dp / (c * c), a3 = (dp + rho * c * du) / (2 * c * c);
            double l1 = std::abs(u - c), l2 = std::abs(u), l3 = std::abs(u + c);
            double r1[3] = {1, u - c, H - u * c}, r2[3] = {1, u, 0.5 * u * u}, r3[3] = {1, u + c, H + u * c};
            for (int i = 0; i < 3; ++i)
                Fn[i][k] = 0.5 * (FL[i] + FR[i]) - 0.5 * (l1 * a1 * r1[i] + l2 * a2 * r2[i] + l3 * a3 * r3[i]);
        }
    }
    void print(const char *name, const Vec &v) const
    {
        if (!trace) return;
        std::printf("%-34s", name);
        for (double x : v) std::printf(" %10.6f", x);
        std::printf("\n");
    }
    // AdvectionWeakDG::Advect -> returns div(F) in weak DG form, in physical space
    Field DoAdvection(const Field &U) const
    {
        int nv = 3, ne = f.ne;
        Field F = GetFluxVector(U);                                       // (1) volume flux
        std::vector<std::vector<Vec>> c(nv, std::vector<Vec>(ne));
        for (int i = 0; i < nv; ++i)
            for (int e = 0; e < ne; ++e) c[i][e] = f.exp.IProductWRTDerivBase(F[i][e]);   // (2)
        std::vector<Vec> Fwd(nv), Bwd(nv), Fn;
        for (int i = 0; i < nv; ++i) f.GetFwdBwdTracePhys(U[i], Fwd[i], Bwd[i]);         // (3)
        RiemannSolve(Fwd, Bwd, Fn);                                                     // (4)
        if (trace)
        {
            for (int e = 0; e < ne; ++e) { print("rho phys, elmt", U[0][e]); }
            for (int e = 0; e < ne; ++e) print("F_rho = rho u at quad pts", F[0][e]);
            for (int e = 0; e < ne; ++e) print("IProductWRTDerivBase(F_rho)", c[0][e]);
            print("Fwd rho (trace)", Fwd[0]); print("Bwd rho (trace)", Bwd[0]); print("Fn rho (Riemann)", Fn[0]);
        }
        for (int i = 0; i < nv; ++i)
            for (int e = 0; e < ne; ++e) for (double &x : c[i][e]) x = -x;                // (5) Neg
        for (int i = 0; i < nv; ++i) f.AddTraceIntegral(Fn[i], c[i]);                   // (6)
        if (trace) for (int e = 0; e < ne; ++e) print("-IProd + AddTraceIntegral", c[0][e]);
        Field out(nv, std::vector<Vec>(ne));
        for (int i = 0; i < nv; ++i)
            for (int e = 0; e < ne; ++e)
                out[i][e] = f.exp.BwdTrans(f.exp.MultiplyByElmtInvMass(c[i][e]));        // (7)(8)
        if (trace) for (int e = 0; e < ne; ++e) print("div F after InvMass + BwdTrans", out[0][e]);
        return out;
    }

    // Orthonormal-Legendre modal filter Q applied element by element (SVV kernel)
    Vec ApplyKernel(const Vec &u) const
    {
        int Q = f.exp.nq, P = f.exp.ncoeffs - 1, M = int(svvCutoff * P);
        Vec out(Q, 0);
        for (int k = 0; k <= P; ++k)
        {
            if (k <= M) continue;
            double kern = std::exp(-double((k - P) * (k - P)) / double((k - M) * (k - M)));
            double uk = 0;
            for (int q = 0; q < Q; ++q) uk += f.exp.w[q] * Jacobi(k, 0, 0, f.exp.z[q]) * u[q];
            uk *= (2 * k + 1) / 2.0;
            for (int q = 0; q < Q; ++q) out[q] += kern * uk * Jacobi(k, 0, 0, f.exp.z[q]);
        }
        return out;
    }
    // weak derivative with a central trace value (LDG/BR1 first step)
    std::vector<Vec> WeakDeriv(const std::vector<Vec> &u, bool zeroWallFlux = false) const
    {
        int ne = f.ne;
        Vec Fwd, Bwd; f.GetFwdBwdTracePhys(u, Fwd, Bwd);
        Vec avg(Fwd.size()); for (size_t k = 0; k < avg.size(); ++k) avg[k] = 0.5 * (Fwd[k] + Bwd[k]);
        if (zeroWallFlux && !f.periodic) { avg[0] = 0; avg[ne] = 0; }   // no SVV flux through the domain ends
        std::vector<Vec> c(ne);
        for (int e = 0; e < ne; ++e) { c[e] = f.exp.IProductWRTDerivBase(u[e]); for (double &x : c[e]) x = -x; }
        f.AddTraceIntegral(avg, c);
        std::vector<Vec> q(ne);
        for (int e = 0; e < ne; ++e) q[e] = f.exp.BwdTrans(f.exp.MultiplyByElmtInvMass(c[e]));
        return q;
    }
    // SVV term d/dx( eps Q[dU/dx] ), eps = SVVDiffCoeff * h/P * lambda_max
    Field DoSVV(const Field &U, double lam) const
    {
        Field out(3, std::vector<Vec>(f.ne));
        double eps = svvCoeff * f.h / (f.exp.ncoeffs - 1) * lam;
        for (int i = 0; i < 3; ++i)
        {
            std::vector<Vec> q = WeakDeriv(U[i]);
            for (int e = 0; e < f.ne; ++e) { q[e] = ApplyKernel(q[e]); for (double &x : q[e]) x *= eps; }
            out[i] = WeakDeriv(q, true);
        }
        return out;
    }
    double MaxSpeed(const Field &U) const
    {
        double lam = 0;
        for (int e = 0; e < f.ne; ++e)
            for (int q = 0; q < f.exp.nq; ++q)
            {
                double r = U[0][e][q], u = U[1][e][q] / r, p = pressure(r, U[1][e][q], U[2][e][q]);
                lam = std::max(lam, std::abs(u) + std::sqrt(GAMMA * std::max(p, 1e-14) / r));
            }
        return lam;
    }
    // CompressibleFlowSystem::DoOdeRhs
    Field DoOdeRhs(const Field &U) const
    {
        Field out = DoAdvection(U);
        for (auto &v : out) for (auto &e : v) for (double &x : e) x = -x;               // Vmath::Neg
        if (svvCoeff > 0)
        {
            Field s = DoSVV(U, MaxSpeed(U));
            for (int i = 0; i < 3; ++i) for (int e = 0; e < f.ne; ++e)
                for (int q = 0; q < f.exp.nq; ++q) out[i][e][q] += s[i][e][q];
        }
        if (trace) for (int e = 0; e < f.ne; ++e) print("d rho/dt = DoOdeRhs", out[0][e]);
        return out;
    }
    // DoOdeProjection: DG keeps physical values as they are (BCs live in the trace here)
    Field DoOdeProjection(const Field &U) const { return U; }

    void RK4(Field &U, double dt) const
    {
        auto axpy = [](const Field &a, const Field &b, double s) {
            Field r = a;
            for (size_t i = 0; i < r.size(); ++i) for (size_t e = 0; e < r[i].size(); ++e)
                for (size_t q = 0; q < r[i][e].size(); ++q) r[i][e][q] += s * b[i][e][q];
            return r; };
        Field k1 = DoOdeRhs(DoOdeProjection(U));
        Field k2 = DoOdeRhs(DoOdeProjection(axpy(U, k1, dt / 2)));
        Field k3 = DoOdeRhs(DoOdeProjection(axpy(U, k2, dt / 2)));
        Field k4 = DoOdeRhs(DoOdeProjection(axpy(U, k3, dt)));
        for (size_t i = 0; i < U.size(); ++i) for (size_t e = 0; e < U[i].size(); ++e)
            for (size_t q = 0; q < U[i][e].size(); ++q)
                U[i][e][q] += dt / 6 * (k1[i][e][q] + 2 * k2[i][e][q] + 2 * k3[i][e][q] + k4[i][e][q]);
    }
    // SetInitialConditions: evaluate, then FwdTrans + BwdTrans (project onto the polynomial space)
    Field SetInitialConditions(std::function<void(double, double &, double &, double &)> ic) const
    {
        Field U(3, std::vector<Vec>(f.ne, Vec(f.exp.nq)));
        for (int e = 0; e < f.ne; ++e)
        {
            for (int q = 0; q < f.exp.nq; ++q)
            {
                double r, u, p; ic(f.xq(e, q), r, u, p);
                U[0][e][q] = r; U[1][e][q] = r * u; U[2][e][q] = p / (GAMMA - 1) + 0.5 * r * u * u;
            }
            for (int i = 0; i < 3; ++i) U[i][e] = f.exp.BwdTrans(f.exp.FwdTrans(U[i][e]));
        }
        return U;
    }
    double GetTimeStep(const Field &U, double cfl) const
    {
        int P = f.exp.ncoeffs - 1;
        double lam = MaxSpeed(U), dt = cfl * f.h / ((2 * P + 1) * lam);
        double eps = svvCoeff * f.h / P * lam;            // explicit viscous limit for the SVV term
        if (eps > 0) dt = std::min(dt, cfl * f.h * f.h / (eps * std::pow(P + 1, 4)));
        return dt;
    }
};

int main(int argc, char **argv)
{
    std::string mode = argc > 1 ? argv[1] : "trace";
    auto wave = [](double x, double &r, double &u, double &p) { r = 1 + 0.2 * std::sin(PI * x); u = 1; p = 1; };

    if (mode == "trace")                          // one DoOdeRhs, 2 elements, P = 2, Q = 4
    {
        EulerCFE s(2, 2, -1, 1, true, "LaxFriedrichs");
        s.trace = true;
        Field U = s.SetInitialConditions(wave);
        Field R = s.DoOdeRhs(U);
        for (int e = 0; e < 2; ++e)
        {
            Vec ex(s.f.exp.nq);
            for (int q = 0; q < s.f.exp.nq; ++q) ex[q] = -0.2 * PI * std::cos(PI * s.f.xq(e, q));
            s.print("exact -d(rho u)/dx", ex);
        }
        Vec xs; for (int e = 0; e < 2; ++e) for (int q = 0; q < s.f.exp.nq; ++q) xs.push_back(s.f.xq(e, q));
        s.print("x (quadrature points)", xs);
        return 0;
    }
    if (mode == "conv")                           // smooth density wave, t = 1 (half period)
    {
        for (std::string up : {"LaxFriedrichs", "Roe"})
            for (int P = 1; P <= 5; ++P)
            {
                std::printf("%-13s P=%d", up.c_str(), P);
                for (int ne : {8, 16})
                {
                    EulerCFE s(ne, P, -1, 1, true, up);
                    Field U = s.SetInitialConditions(wave);
                    double T = 1.0, dt = s.GetTimeStep(U, 0.2);
                    int n = std::ceil(T / dt); dt = T / n;
                    for (int k = 0; k < n; ++k) s.RK4(U, dt);
                    double err = 0;
                    for (int e = 0; e < ne; ++e) for (int q = 0; q < s.f.exp.nq; ++q)
                        err = std::max(err, std::abs(U[0][e][q] - (1 + 0.2 * std::sin(PI * (s.f.xq(e, q) - T)))));
                    std::printf("  Ne=%2d err=%.3e", ne, err);
                }
                std::printf("\n");
            }
        return 0;
    }
    if (mode == "profile")                       // rho(x) at t = 0.2 for SVVDiffCoeff = 0.5
    {
        auto sod = [](double x, double &r, double &u, double &p) {
            if (x < 0.5) { r = 1; u = 0; p = 1; } else { r = 0.125; u = 0; p = 0.1; } };
        EulerCFE s(40, 4, 0, 1, false, "Roe");
        s.svvCoeff = 0.5;
        Field U = s.SetInitialConditions(sod);
        double t = 0;
        while (t < 0.2 - 1e-12) { double dt = std::min(s.GetTimeStep(U, 0.1), 0.2 - t); s.RK4(U, dt); t += dt; }
        for (int e = 0; e < 40; ++e) for (int q = 0; q < s.f.exp.nq; ++q) std::printf("%.5f %.5f\n", s.f.xq(e, q), U[0][e][q]);
        return 0;
    }
    if (mode == "sod")                            // Sod tube, P = 4, 40 elements, t = 0.2
    {
        auto sod = [](double x, double &r, double &u, double &p) {
            if (x < 0.5) { r = 1; u = 0; p = 1; } else { r = 0.125; u = 0; p = 0.1; } };
        for (double c : {0.0, 0.25, 0.5, 1.0})
        {
            EulerCFE s(40, 4, 0, 1, false, "Roe");
            s.svvCoeff = c;
            Field U = s.SetInitialConditions(sod);
            double t = 0; bool ok = true;
            while (t < 0.2 - 1e-12)
            {
                double dt = std::min(s.GetTimeStep(U, 0.1), 0.2 - t);
                s.RK4(U, dt); t += dt;
                for (auto &el : U[0]) for (double r : el) if (!std::isfinite(r) || r <= 0) ok = false;
                if (!ok) break;
            }
            double rmax = 0, rmin = 1e9, tv = 0, prev = U[0][0][0];
            for (int e = 0; e < 40; ++e) for (int q = 0; q < s.f.exp.nq; ++q)
            {
                double r = U[0][e][q];
                rmax = std::max(rmax, r); rmin = std::min(rmin, r); tv += std::abs(r - prev); prev = r;
            }
            if (ok) std::printf("SVVDiffCoeff=%.2f  t=0.2  max rho=%.4f  min rho=%.4f  TV(rho)=%.4f\n", c, rmax, rmin, tv);
            else    std::printf("SVVDiffCoeff=%.2f  failed (non-finite state) at t=%.4f\n", c, t);
        }
        return 0;
    }
}
