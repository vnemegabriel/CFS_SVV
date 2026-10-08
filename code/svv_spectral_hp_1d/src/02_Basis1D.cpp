#include "02_Basis1D.h"
#include "polylib.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

// -----------------------------------------------------------------------------
// Modified basis (K&S 2.40).  polylib::jacobfd gives P^{a,b}_n(z), and
// polylib::jacobd gives its derivative in a form valid at z = +-1 as well.
// -----------------------------------------------------------------------------
void Basis1D::evalModifiedBasis(int P, const Vector &z, Matrix &B, Matrix &dB)
{
    const int np = static_cast<int>(z.size());
    B            = Matrix(np, P + 1);
    dB           = Matrix(np, P + 1);

    Vector zz(z);  // polylib wants non-const pointers
    Vector J(np), dJ(np);

    for (int q = 0; q < np; ++q)
    {
        B(q, 0)  = 0.5 * (1.0 - z[q]);
        dB(q, 0) = -0.5;
        B(q, P)  = 0.5 * (1.0 + z[q]);
        dB(q, P) = 0.5;
    }
    for (int p = 1; p <= P - 1; ++p)
    {
        polylib::jacobfd(np, zz.data(), J.data(), nullptr, p - 1, 1.0, 1.0);
        polylib::jacobd(np, zz.data(), dJ.data(), p - 1, 1.0, 1.0);
        for (int q = 0; q < np; ++q)
        {
            const double omz2 = 1.0 - z[q] * z[q];
            B(q, p)           = 0.25 * omz2 * J[q];
            //  d/dxi [ (1-xi^2)/4 * J ] = -xi/2 * J + (1-xi^2)/4 * J'
            dB(q, p) = -0.5 * z[q] * J[q] + 0.25 * omz2 * dJ[q];
        }
    }
}

// -----------------------------------------------------------------------------
// Orthonormal Legendre:  psi~_k = sqrt((2k+1)/2) L_k,  L_k = P^{0,0}_k
// -----------------------------------------------------------------------------
Matrix Basis1D::evalOrthonormalLegendre(int P, const Vector &z)
{
    const int np = static_cast<int>(z.size());
    Matrix L(np, P + 1);
    Vector zz(z), Lk(np);
    for (int k = 0; k <= P; ++k)
    {
        polylib::jacobfd(np, zz.data(), Lk.data(), nullptr, k, 0.0, 0.0);
        const double nrm = std::sqrt((2.0 * k + 1.0) / 2.0);
        for (int q = 0; q < np; ++q)
            L(q, k) = nrm * Lk[q];
    }
    return L;
}

// -----------------------------------------------------------------------------
Basis1D::Basis1D(int P_, int Q_) : P(P_), nModes(P_ + 1), Q(Q_)
{
    if (P < 1)
        throw std::runtime_error("Basis1D: need P >= 1");
    if (Q < P + 1)
        throw std::runtime_error("Basis1D: need Q >= P+1 quadrature points");

    // --- solution quadrature -------------------------------------------------
    zq.assign(Q, 0.0);
    wq.assign(Q, 0.0);
    polylib::zwglj(zq.data(), wq.data(), Q, 0.0, 0.0);
    evalModifiedBasis(P, zq, B, dB);
    Bortho = evalOrthonormalLegendre(P, zq);

    // --- exact quadrature for the reference matrices ---------------------------
    // GLL with Qe points is exact to degree 2Qe-3; we need degree 2P => Qe = P+2.
    const int Qe = P + 2;
    Vector ze(Qe), we(Qe);
    polylib::zwglj(ze.data(), we.data(), Qe, 0.0, 0.0);
    Matrix Be, dBe;
    evalModifiedBasis(P, ze, Be, dBe);
    Matrix Boe = evalOrthonormalLegendre(P, ze);

    // (f, g) = sum_q w_q f(z_q) g(z_q)   ->   F^T W G
    const Matrix We = Matrix::diagonal(we);
    M                = Be.transpose() * We * Be;
    S                = Be.transpose() * We * dBe;
    K                = dBe.transpose() * We * dBe;
    T                = Boe.transpose() * We * Be;
    Minv             = M.inverse();
    Tinv             = T.inverse();
}

Vector Basis1D::innerProduct(const Vector &fq) const
{
    return B.transpose() * hadamard(wq, fq);
}

Vector Basis1D::innerProductDeriv(const Vector &fq) const
{
    return dB.transpose() * hadamard(wq, fq);
}

double Basis1D::minGLLSpacing() const
{
    Vector z(P + 1), w(P + 1);
    polylib::zwglj(z.data(), w.data(), P + 1, 0.0, 0.0);
    double dmin = 2.0;
    for (int i = 1; i <= P; ++i)
        dmin = std::min(dmin, z[i] - z[i - 1]);
    return dmin;
}
