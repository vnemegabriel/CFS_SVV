// =============================================================================
// 08_Burgers.h  --  inviscid Burgers + SVV, the K&S Fig. 6.27 test
//
//      u_t + (u^2/2)_x = eps d/dx [ F * u_x ]               (6.5.14)
//
// Weak form with a C0 test function v (integrating the SVV term by parts):
//
//      (v, u_t) = -(v, (u^2/2)_x) - eps (v_x, F u_x)
//
// Discretisation per element (reference coordinates, J = h/2):
//      mass       :  J M_ref
//      convection :  (psi_i, df/dx) = S_ref f^ ,    f^ = M^{-1} B^T W (u_q^2/2)
//                    (the 2/h from d/dx cancels the J of the integral)
//      SVV        :  (eps/J) L_svv,ref u^
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"
#include "03_SVV.h"
#include "05_CGAssembly.h"

enum class ConvectiveForm
{
    Conservative,      // (v, (u^2/2)_x)
    NonConservative,   // (v, u u_x)
    Skew               // (1/3)(v, (u^2)_x) + (1/3)(v, u u_x)
};


struct BurgersOptions
{
    ConvectiveForm conv = ConvectiveForm::Conservative;
    SVVApplication apply = SVVApplication::Split;
    bool useSVV      = true;
    double eps       = 1.0 / 16.0;   // K&S: eps = 1/16
    int Mcut         = 8;            // K&S: M_SVV = 8
    SVVKernelType kernel = SVVKernelType::Exponential;
    SVVForm form     = SVVForm::Book;
};

class BurgersSolver
{
public:
    BurgersSolver(const CGAssembly &A, const BurgersOptions &opt);

    // R(u^) = M_bc^{-1} r(u^) -- what the time integrator sees
    Vector rhs(const Vector &uhat, double t) const;

    // Split application: call once after every time step of size dt.
    void applySplitSVV(Vector &uhat, double dt);

    // Largest |eigenvalue| of M^{-1} (eps L_svv): sets the diffusive dt limit
    double svvSpectralRadius() const;

    const Matrix &svvReference() const { return Lref_; }
    const Vector &kernel() const { return kernel_; }

private:
    const CGAssembly &A_;
    BurgersOptions opt_;
    Vector kernel_;
    Matrix Lref_;      // reference SVV operator
    LUSolver massLU_;  // mass matrix with Dirichlet rows replaced
    SVVBubbleFilter bubbleFilter_;   // split application (declared last: initialised last)
};
