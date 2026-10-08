#include "08_Burgers.h"

namespace
{
Matrix buildMassBC(const CGAssembly &A)
{
    Matrix M = A.massMatrix();
    A.applyDirichletRows(M);
    return M;
}
}   // namespace

BurgersSolver::BurgersSolver(const CGAssembly &A, const BurgersOptions &opt)
    : A_(A), opt_(opt), kernel_(svvKernel(A.P, opt.Mcut, opt.kernel)),
      Lref_(svvOperator(A.basis, kernel_, opt.form)), massLU_(buildMassBC(A)), bubbleFilter_(A.basis, Lref_)
{
}

void BurgersSolver::applySplitSVV(Vector &uhat, double dt)
{
    if (!opt_.useSVV || opt_.apply != SVVApplication::Split)
        return;
    for (int e = 0; e < A_.nel; ++e)
    {
        // element system (J M) du/dt = -(eps/J) L u, divided by J -> c = dt eps / J^2
        const double J  = A_.mesh.J(e);
        const Vector ue = bubbleFilter_.apply(A_.gather(uhat, e), dt * opt_.eps / (J * J));
        for (int p = 1; p < A_.P; ++p)   // bubble dofs are not shared: plain overwrite
            uhat[A_.dof(e, p)] = ue[p];
    }
}

Vector BurgersSolver::rhs(const Vector &uhat, double) const
{
    const Basis1D &b = A_.basis;
    Vector r(A_.nGlobal, 0.0);
    Vector fq(b.Q);

    for (int e = 0; e < A_.nel; ++e)
    {
        const Vector ue = A_.gather(uhat, e);

        // --- convection ---------------------------------------------------------
        const Vector uq = b.evaluate(ue);
        Vector re;
        if (opt_.conv == ConvectiveForm::Conservative)
        {
            // (psi_i, d/dx (u^2/2)) = S f^,  f^ = projection of u^2/2 (aliased unless Q is large)
            for (int q = 0; q < b.Q; ++q)
                fq[q] = 0.5 * uq[q] * uq[q];
            re = scaled(b.S * b.project(fq), -1.0);
        }
        else if (opt_.conv == ConvectiveForm::NonConservative)
        {
            // (psi_i, u du/dx) = B^T W (u .* dB u^)   (J and 1/J cancel)
            re = scaled(b.innerProduct(hadamard(uq, b.evaluateDeriv(ue))), -1.0);
        }
        else
        {
            // (1/3)(psi_i, d/dx u^2) + (1/3)(psi_i, u du/dx)
            // second term: J * (psi_i, u (1/J) du/dxi)_ref = B^T W (u .* dB u^)
            for (int q = 0; q < b.Q; ++q)
                fq[q] = uq[q] * uq[q];
            re                = scaled(b.S * b.project(fq), -1.0 / 3.0);
            const Vector udux = hadamard(uq, b.evaluateDeriv(ue));
            axpy(-1.0 / 3.0, b.innerProduct(udux), re);
        }

        // --- SVV (Galerkin): -(eps/J) L_ref u^ as a weak-form residual ----------
        if (opt_.useSVV && opt_.apply == SVVApplication::Galerkin)
            axpy(-opt_.eps / A_.mesh.J(e), Lref_ * ue, re);

        A_.scatterAdd(re, e, r);
    }

    A_.zeroDirichlet(r);      // du/dt = 0 on Dirichlet dofs
    return massLU_.solve(r);
}

double BurgersSolver::svvSpectralRadius() const
{
    if (!opt_.useSVV)
        return 0.0;
    return opt_.eps * A_.spectralRadius(A_.operatorMatrix(Lref_));
}
