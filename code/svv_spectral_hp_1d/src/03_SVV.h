// =============================================================================
// 03_SVV.h  --  Spectral Vanishing Viscosity for the continuous Galerkin
//               spectral/hp method, K&S section 6.5.2 (pp. 348-352)
//
// CONTINUOUS IDEA (Tadmor 1989, eq. 6.5.10 / 6.5.11)
//     u_t + f(u)_x = eps d/dx [ F * du/dx ]
// The convolution kernel F acts only on the high modes: in Fourier space it is
// a multiplication of the coefficients of u_x by F^_k, with F^_k = 0 for
// |k| <= M_SVV.  Low (well resolved) modes therefore see NO viscosity and the
// method keeps spectral accuracy; the high modes see an O(eps) viscosity that
// kills the wiggles.
//
// KERNEL (Maday, Ould Kaber & Tadmor 1993, eq. 6.5.13)
//     F^_k = exp( -(k-P)^2 / (k-M_SVV)^2 )   for k > M_SVV,   0 otherwise.
//     -> F^_P = 1 (full viscosity on the last mode), C-infinity ramp in between.
//
// DISCRETE OPERATOR FOR C0 GALERKIN (Kirby 2002, K&S p. 351)
// The weak form of the SVV term is  ( dv/dx, F du/dx ).  Per element:
//     1. u^        -> modal coefficients of u^' = du/dxi:      w^  = M^{-1} S u^
//     2. filter in the ORTHONORMAL space (where "mode = degree"):
//                                                           w^f = T^{-1} F T w^
//     3. inner product with dv/dxi:                          S^T w^f
// so that
//     SVV_ref = S^T T^{-1} F T M^{-1} S  =  S^T M^{-1} T^T F T M^{-1} S     (*)
// symmetric, positive semi-definite. Physical scaling: (dv/dx, F du/dx) on an
// element of length h is (2/h) * SVV_ref  (one 2/h per derivative, times the
// Jacobian h/2).  With F = I the operator reduces to the Laplacian K.
//
// ALTERNATIVE (Kirby, eq. 6.5.17): SVV2 = eps F d2/dx2 F u, i.e. filter the
// solution, take the Laplacian, filter again:  Phi^T K Phi,  Phi = T^{-1} F T.
// Less dissipative than (*).
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"
#include "02_Basis1D.h"

#include <memory>
#include <string>

enum class SVVKernelType
{
    Exponential,   // eq. 6.5.13 (default)
    Step           // Tadmor's original 0/1 kernel, eq. 6.5.12
};

enum class SVVForm
{
    Book,    // S^T T^{-1} F T M^{-1} S           (eq. on p. 351)
    Kirby    // Phi^T K Phi, Phi = T^{-1} F T     (eq. 6.5.17)
};

SVVKernelType parseKernelType(const std::string &s);
SVVForm parseSVVForm(const std::string &s);

// F^_k for k = 0..P  (diagonal of the filter matrix F).
Vector svvKernel(int P, int Mcut, SVVKernelType type = SVVKernelType::Exponential);

// Reference-element SVV matrix (nModes x nModes). Multiply by (2/h) = 1/J for
// an element of size h.
Matrix svvOperator(const Basis1D &basis, const Vector &kernel, SVVForm form = SVVForm::Book);

// The modal filter Phi = T^{-1} F T acting on modified-basis coefficients.
Matrix modalFilter(const Basis1D &basis, const Vector &kernel);

// ----------------------------------------------------------------------------
// TWO WAYS OF APPLYING THE SVV OPERATOR IN A C0 (continuous Galerkin) CODE
//
//  Galerkin (the formula on p. 351, used as a weak-form residual):
//        M_global du^/dt = r_conv(u^) - sum_e (eps/J_e) L_e u^_e
//     The vertex rows/columns of L_e are ZERO (the derivative of a vertex mode
//     is constant, i.e. Legendre mode 0, which every kernel annihilates), so
//     SVV produces residuals on bubble modes only. The consistent mass matrix
//     nevertheless couples those residuals to the shared vertex dof: element
//     e's SVV moves the vertex, and the neighbour must absorb a moving vertex
//     with a function L2-orthogonal to its bubbles -- a spiky, high-degree
//     "lift". In an inviscid problem nothing damps that lift; the result is a
//     spike at every element interface (see docs/03_burgers_test.md).
//
//  Split (element-wise implicit filter, vertex frozen):
//     after every explicit step, solve in each element the bubble subsystem of
//        (J_e M) (u^{n+1} - u^n) = -dt (eps/J_e) L u^{n+1}, divided by J_e:
//        (M_bb + dt (eps/J_e^2) L_bb) u_b^{n+1} = M_bb u_b^n
//     (b = bubble modes 1..P-1). Because L has no vertex coupling this is the
//     exact bubble block of the implicit-Euler SVV step with the vertex held
//     fixed. It is unconditionally stable, keeps C0 continuity trivially, never
//     talks to the neighbours, and is what a spectral/hp code does when it
//     "filters the interior modes".
// ----------------------------------------------------------------------------
enum class SVVApplication
{
    Galerkin,
    Split
};
SVVApplication parseSVVApplication(const std::string &s);

// Element-wise implicit bubble filter, precomputed for a given coefficient
// c = dt * eps / J^2 (reference mass on the left, hence J^2).
// Rebuilds its LU only when c changes.
class SVVBubbleFilter
{
public:
    SVVBubbleFilter(const Basis1D &basis, const Matrix &Lref);
    // u^_e -> filtered u^_e (vertex coefficients returned unchanged)
    Vector apply(const Vector &ue, double c);

private:
    const Basis1D &basis_;
    Matrix Mbb_, Lbb_;
    double c_ = -1.0;
    std::unique_ptr<LUSolver> lu_;
};
