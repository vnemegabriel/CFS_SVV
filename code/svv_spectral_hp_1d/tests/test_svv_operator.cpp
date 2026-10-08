// =============================================================================
// Unit tests: every identity quoted from the book is checked numerically.
// Run:  ./build/test_svv_operator      (or: ctest --test-dir build)
// =============================================================================
#include "02_Basis1D.h"
#include "03_SVV.h"
#include "04_Mesh1D.h"
#include "05_CGAssembly.h"
#include "08_Burgers.h"

#include <cmath>
#include <cstdio>
#include <functional>

static int nPass = 0, nFail = 0;
#define CHECK(cond, msg)                                                                                     \
    do                                                                                                       \
    {                                                                                                        \
        if (cond)                                                                                            \
        {                                                                                                    \
            ++nPass;                                                                                         \
            std::printf("  PASS  %s\n", msg);                                                                \
        }                                                                                                    \
        else                                                                                                 \
        {                                                                                                    \
            ++nFail;                                                                                         \
            std::printf("  FAIL  %s\n", msg);                                                                \
        }                                                                                                    \
    } while (0)

// modal coefficients of the Legendre polynomial of degree j (exact projection)
static Vector legendreModeCoeffs(const Basis1D &b, int j)
{
    Vector fq(b.Q);
    for (int q = 0; q < b.Q; ++q)
        fq[q] = b.Bortho(q, j);
    return b.project(fq);
}

int main()
{
    const int P = 15, Mcut = 8;
    Basis1D b(P, P + 2);
    const double tol = 1e-9;

    std::printf("[basis identities, P=%d]\n", P);
    CHECK((b.T.transpose() * b.T - b.M).maxAbs() < tol, "Parseval: T^T T = M");
    CHECK((b.Tinv - b.Minv * b.T.transpose()).maxAbs() < tol, "T^{-1} = M^{-1} T^T (K&S p. 351)");
    CHECK((b.S.transpose() * b.Minv * b.S - b.K).maxAbs() < tol, "S^T M^{-1} S = K (derivative projection is exact)");
    CHECK(b.M.isSymmetric(tol) && b.K.isSymmetric(tol), "M and K symmetric");
    {
        // vertex modes are the only ones non-zero at xi = +-1
        Matrix B, dB;
        Basis1D::evalModifiedBasis(P, {-1.0, 1.0}, B, dB);
        bool ok = std::fabs(B(0, 0) - 1.0) < tol && std::fabs(B(1, P) - 1.0) < tol;
        for (int p = 1; p < P; ++p)
            ok = ok && std::fabs(B(0, p)) < tol && std::fabs(B(1, p)) < tol;
        CHECK(ok, "only vertex modes are non-zero on the element boundary");
    }
    {
        // derivative check: project x^5 (exact for P >= 5) and differentiate -> 5 x^4
        Vector fq(b.Q);
        for (int q = 0; q < b.Q; ++q)
            fq[q] = std::pow(b.zq[q], 5);
        const Vector du = b.evaluateDeriv(b.project(fq));
        double err      = 0.0;
        for (int q = 0; q < b.Q; ++q)
            err = std::max(err, std::fabs(du[q] - 5.0 * std::pow(b.zq[q], 4)));
        CHECK(err < 1e-9, "dB differentiates the modal expansion correctly (d/dxi x^5 = 5 x^4)");
    }

    std::printf("[SVV kernel]\n");
    const Vector F = svvKernel(P, Mcut);
    {
        bool zeroLow = true, mono = true;
        for (int k = 0; k <= Mcut; ++k)
            zeroLow = zeroLow && F[k] == 0.0;
        for (int k = Mcut + 1; k < P; ++k)
            mono = mono && F[k + 1] > F[k];
        CHECK(zeroLow, "F_k = 0 for k <= M_SVV");
        CHECK(std::fabs(F[P] - 1.0) < 1e-15, "F_P = 1");
        CHECK(mono, "F_k increases monotonically for k > M_SVV");
    }

    std::printf("[SVV operator]\n");
    const Matrix L  = svvOperator(b, F, SVVForm::Book);
    const Matrix L2 = svvOperator(b, F, SVVForm::Kirby);
    CHECK(L.isSymmetric(tol), "book form is symmetric");
    CHECK(L2.isSymmetric(tol), "Kirby form is symmetric");
    {
        bool psd = true, psd2 = true;
        for (int trial = 0; trial < 20; ++trial)
        {
            Vector v(b.nModes);
            for (int i = 0; i < b.nModes; ++i)
                v[i] = std::sin(3.1 * i + trial);
            psd  = psd && dot(v, L * v) >= -tol;
            psd2 = psd2 && dot(v, L2 * v) >= -tol;
        }
        CHECK(psd, "book form is positive semi-definite");
        CHECK(psd2, "Kirby form is positive semi-definite");
    }
    CHECK((svvOperator(b, Vector(b.nModes, 1.0)) - b.K).maxAbs() < tol, "F = I recovers the Laplacian K");
    {
        const Matrix MPhi = b.M * modalFilter(b, F);
        CHECK(MPhi.isSymmetric(tol), "modal filter is L2-symmetric: M Phi = T^T F T");
        bool zeroVertex = true;
        for (int j = 0; j < b.nModes; ++j)
            for (int v : {0, P})
                zeroVertex = zeroVertex && std::fabs(L(v, j)) < tol && std::fabs(L(j, v)) < tol;
        CHECK(zeroVertex, "vertex rows and columns of L_svv are zero");
    }
    {
        // SVV acts on u_x. If deg(u) <= M+1 then deg(u_x) <= M and the filter kills it.
        bool vanish = true;
        for (int j = 0; j <= Mcut + 1; ++j)
            vanish = vanish && maxAbs(L * legendreModeCoeffs(b, j)) < 1e-8;
        CHECK(vanish, "book form annihilates polynomials of degree <= M_SVV+1 (no viscosity on resolved modes)");
        CHECK(maxAbs(L * legendreModeCoeffs(b, P)) > 1.0, "book form acts on the highest mode");
        bool vanish2 = true;
        for (int j = 0; j <= Mcut; ++j)
            vanish2 = vanish2 && maxAbs(L2 * legendreModeCoeffs(b, j)) < 1e-8;
        CHECK(vanish2, "Kirby form annihilates polynomials of degree <= M_SVV");
    }
    {
        // Energy argument: with F = I the SVV term dissipates exactly like a Laplacian;
        // with the kernel it dissipates less:  u^T L u <= u^T K u.
        Vector v(b.nModes);
        for (int i = 0; i < b.nModes; ++i)
            v[i] = std::cos(1.3 * i);
        CHECK(dot(v, L * v) <= dot(v, b.K * v) + tol, "SVV dissipation <= Laplacian dissipation");
    }

    std::printf("[split application: element-wise implicit bubble filter]\n");
    {
        SVVBubbleFilter filt(b, L);
        const double c = 0.05;   // dt * eps / J^2
        // low-degree content passes untouched, vertex coefficients never change
        bool passLow = true;
        for (int j = 0; j <= Mcut + 1; ++j)
        {
            const Vector u  = legendreModeCoeffs(b, j);
            const Vector uf = filt.apply(u, c);
            for (int p = 0; p < b.nModes; ++p)
                passLow = passLow && std::fabs(uf[p] - u[p]) < 1e-8;
        }
        CHECK(passLow, "bubble filter leaves polynomials of degree <= M_SVV+1 unchanged");
        Vector top       = legendreModeCoeffs(b, P);
        const Vector tf  = filt.apply(top, c);
        const Vector ttf = b.toOrthonormal(tf);
        CHECK(std::fabs(ttf[P]) < 0.5 * std::fabs(b.toOrthonormal(top)[P]), "bubble filter damps the highest mode");
        CHECK(tf[0] == top[0] && tf[P] == top[P], "bubble filter never touches the vertex coefficients");
        // energy never increases: u^T M u_f <= u^T M u
        CHECK(dot(tf, b.M * tf) <= dot(top, b.M * top) + tol, "bubble filter is dissipative in the L2 norm");
    }

    std::printf("[split application: physical scaling]\n");
    {
        // One split step must equal the bubble block of implicit Euler on the PHYSICAL
        // element system (J M)(u1 - u0) = -dt (eps/J) L u1 with the vertices frozen.
        Mesh1D mesh = Mesh1D::uniform(-1.0, 1.0, 5);   // J = 0.2
        CGAssembly A(mesh, b, BCType::Dirichlet);
        BurgersOptions opt;
        opt.apply = SVVApplication::Split;
        BurgersSolver solver(A, opt);
        Vector u(A.nGlobal);
        for (int i = 0; i < A.nGlobal; ++i)
            u[i] = std::cos(1.7 * i);
        Vector us       = u;
        const double dt = 2.5e-3;
        solver.applySplitSVV(us, dt);
        const Matrix &Lr = solver.svvReference();
        double err       = 0.0;
        for (int e = 0; e < A.nel; ++e)
        {
            const double J = mesh.J(e);
            const Matrix Mbb = b.M.block(1, P, 1, P).scaled(J);
            const Matrix Lbb = Lr.block(1, P, 1, P).scaled(dt * opt.eps / J);
            const Vector ue  = A.gather(u, e);
            const Vector ub(ue.begin() + 1, ue.begin() + P);
            const Vector ex  = LUSolver(Mbb + Lbb).solve(Mbb * ub);
            for (int p = 1; p < P; ++p)
                err = std::max(err, std::fabs(us[A.dof(e, p)] - ex[p - 1]));
        }
        CHECK(err < 1e-12, "split SVV step = implicit Euler of the physical element system (J scaling)");
    }

    std::printf("[assembly]\n");
    {
        Mesh1D mesh = Mesh1D::uniform(-1.0, 1.0, 5);
        CGAssembly A(mesh, b, BCType::Dirichlet);
        CHECK(A.nGlobal == 5 * P + 1, "Dirichlet: nGlobal = nel*P + 1");
        CHECK(A.dof(0, P) == A.dof(1, 0), "neighbouring elements share the vertex dof");
        // u == 1 has coefficient 1 on every vertex mode and 0 on the bubbles
        // (the modal basis is NOT a partition of unity, so do not sum M blindly).
        Matrix M = A.massMatrix();
        Vector one(A.nGlobal, 0.0);
        for (int e = 0; e < A.nel; ++e)
            one[A.dof(e, 0)] = one[A.dof(e, P)] = 1.0;
        CHECK(std::fabs(dot(one, M * one) - 2.0) < tol, "1^T M 1 = domain length (u == 1 is represented exactly)");
        std::function<double(double)> cubic = [](double x) { return x * x * x - 0.3 * x; };
        Vector uh = A.l2Projection(cubic);
        CHECK(A.l2Norm(uh, &cubic) < 1e-10, "L2 projection of a polynomial is exact");

        CGAssembly Ap(mesh, b, BCType::Periodic);
        CHECK(Ap.nGlobal == 5 * P && Ap.dof(4, P) == 0, "periodic: last vertex is dof 0");

        // spectral radius of M^{-1} K must be positive and O(P^4/h^2)
        const double rho = A.spectralRadius(A.laplacianMatrix());
        CHECK(rho > 1e3 && rho < 1e7, "spectral radius of M^{-1}K is in the expected range");
    }

    std::printf("\n%d passed, %d failed\n", nPass, nFail);
    return nFail == 0 ? 0 : 1;
}
