#include "03_SVV.h"

#include <cmath>
#include <stdexcept>

SVVKernelType parseKernelType(const std::string &s)
{
    if (s == "exp" || s == "exponential")
        return SVVKernelType::Exponential;
    if (s == "step" || s == "tadmor")
        return SVVKernelType::Step;
    throw std::runtime_error("unknown SVV kernel: " + s + " (use exp|step)");
}

SVVForm parseSVVForm(const std::string &s)
{
    if (s == "book" || s == "1")
        return SVVForm::Book;
    if (s == "kirby" || s == "2")
        return SVVForm::Kirby;
    throw std::runtime_error("unknown SVV form: " + s + " (use book|kirby)");
}

SVVApplication parseSVVApplication(const std::string &s)
{
    if (s == "galerkin")
        return SVVApplication::Galerkin;
    if (s == "split")
        return SVVApplication::Split;
    throw std::runtime_error("unknown SVV application: " + s + " (use galerkin|split)");
}

SVVBubbleFilter::SVVBubbleFilter(const Basis1D &basis, const Matrix &Lref) : basis_(basis)
{
    const int P = basis.P;
    Mbb_        = basis.M.block(1, P, 1, P);
    Lbb_        = Lref.block(1, P, 1, P);
}

Vector SVVBubbleFilter::apply(const Vector &ue, double c)
{
    const int P = basis_.P;
    if (P < 2)
        return ue;   // no bubble modes
    if (c != c_)
    {
        lu_ = std::make_unique<LUSolver>(Mbb_ + Lbb_.scaled(c));
        c_  = c;
    }
    Vector ub(ue.begin() + 1, ue.begin() + P);
    const Vector ubNew = lu_->solve(Mbb_ * ub);
    Vector out(ue);
    for (int p = 1; p < P; ++p)
        out[p] = ubNew[p - 1];
    return out;
}

Vector svvKernel(int P, int Mcut, SVVKernelType type)
{
    Vector F(P + 1, 0.0);
    for (int k = Mcut + 1; k <= P; ++k)
    {
        switch (type)
        {
            case SVVKernelType::Exponential:
            {
                // eq. (6.5.13):  exp(-(k-P)^2 / (k-M)^2)
                const double num = static_cast<double>(k - P);
                const double den = static_cast<double>(k - Mcut);
                F[k]             = std::exp(-(num * num) / (den * den));
                break;
            }
            case SVVKernelType::Step:
                F[k] = 1.0;   // eq. (6.5.12)
                break;
        }
    }
    return F;
}

Matrix modalFilter(const Basis1D &basis, const Vector &kernel)
{
    // Phi = T^{-1} F T : modified coeffs -> orthonormal -> scale -> back
    const Matrix F = Matrix::diagonal(kernel);
    return basis.Tinv * F * basis.T;
}

Matrix svvOperator(const Basis1D &basis, const Vector &kernel, SVVForm form)
{
    if (static_cast<int>(kernel.size()) != basis.nModes)
        throw std::runtime_error("svvOperator: kernel length must be P+1");

    const Matrix Phi = modalFilter(basis, kernel);

    switch (form)
    {
        case SVVForm::Book:
        {
            // S^T  Phi  M^{-1} S     with Phi = T^{-1} F T
            //  ^     ^     ^
            //  |     |     +-- coefficients of du/dxi
            //  |     +-------- filter those coefficients
            //  +-------------- inner product with dv/dxi
            const Matrix St = basis.S.transpose();
            return St * Phi * basis.Minv * basis.S;
        }
        case SVVForm::Kirby:
        {
            // (d(Phi v)/dxi, d(Phi u)/dxi) = Phi^T K Phi
            return Phi.transpose() * basis.K * Phi;
        }
    }
    throw std::runtime_error("svvOperator: unreachable");
}
