// =============================================================================
// 02_Basis1D.h  --  the 1D spectral/hp element on the reference interval [-1,1]
//
// Everything in this file is Chapter 2 of Karniadakis & Sherwin (K&S).
//
// TWO BASES FOR THE SAME SPACE  (this is the heart of the SVV construction)
//
//  (1) The modified C0 "p-type" modal basis, K&S eq. (2.40):
//
//        psi_0(xi) = (1 - xi)/2                       left  vertex mode
//        psi_P(xi) = (1 + xi)/2                       right vertex mode
//        psi_p(xi) = (1-xi)/2 (1+xi)/2 P^{1,1}_{p-1}(xi),  1 <= p <= P-1   bubble modes
//
//      Vertex modes are the only ones that are non-zero on the element boundary,
//      so C0 continuity between elements is just "share the vertex coefficient".
//      This basis is NOT orthogonal (semi-orthogonal, K&S p. 351).
//
//  (2) The orthonormal Legendre basis
//
//        psi~_k(xi) = sqrt((2k+1)/2) L_k(xi),   k = 0..P,   (psi~_k, psi~_l) = delta_kl
//
//      Here "mode k" really means "polynomial degree k", which is the notion of
//      frequency that the SVV kernel needs. In Nektar++ this is `eOrtho_A`.
//
// Both span the polynomials of degree <= P on [-1,1]. The change of basis is
//
//        u~ = T u^ ,     T_kp = (psi~_k, psi_p)         (K&S p. 351, "T")
//
// because u~_k = (psi~_k, u) = sum_p u^_p (psi~_k, psi_p).  Note that T^T T = M
// (Parseval), hence T^{-1} = M^{-1} T^T -- the identity used on p. 351.
//
// QUADRATURE
//   * All reference matrices (M, S, K, T) are integrated EXACTLY with a private
//     Gauss-Lobatto-Legendre rule of P+2 points (exact to degree 2P+1).
//   * Non-linear terms (u^2, fluxes) are evaluated on a user-chosen GLL rule of
//     Q >= P+1 points. Q = P+1 is the classic choice; Q > P+1 is "over-
//     integration"/de-aliasing (K&S section 6.5.1, Nektar++ SPECTRALHPDEALIASING).
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"

class Basis1D
{
public:
    // P: polynomial order (number of modes = P+1). Q: quadrature points (>= P+1).
    Basis1D(int P, int Q);

    int P;       // polynomial order
    int nModes;  // P + 1
    int Q;       // number of solution-quadrature points

    // Solution quadrature (GLL) -------------------------------------------
    Vector zq, wq;   // points and weights on [-1,1]
    Matrix B;        // Q x nModes,  B(q,p)  = psi_p(zq[q])
    Matrix dB;       // Q x nModes, dB(q,p)  = d psi_p/d xi (zq[q])
    Matrix Bortho;   // Q x nModes, orthonormal Legendre at the quadrature points

    // Reference-element matrices (exact integration) ------------------------
    Matrix M;        // mass       M_ij = (psi_i, psi_j)
    Matrix Minv;
    Matrix S;        // "stiffness" of K&S p. 351:  S_ij = (psi_i, d psi_j/d xi)
    Matrix K;        // Laplacian  K_ij = (d psi_i/d xi, d psi_j/d xi)
    Matrix T;        // T_kp = (psi~_k, psi_p)   modal -> orthonormal
    Matrix Tinv;     // (equals Minv * T^T because integration is exact)

    // Element-level operations in reference coordinates --------------------
    Vector evaluate(const Vector &uhat) const { return B * uhat; }        // u(zq)
    Vector evaluateDeriv(const Vector &uhat) const { return dB * uhat; }  // du/dxi(zq)
    Vector innerProduct(const Vector &fq) const;        // (psi_i, f)   via quadrature
    Vector innerProductDeriv(const Vector &fq) const;   // (psi_i', f)  via quadrature
    Vector project(const Vector &fq) const { return Minv * innerProduct(fq); }
    Vector toOrthonormal(const Vector &uhat) const { return T * uhat; }

    // Smallest GLL spacing on the P+1 grid: the length scale that governs the
    // explicit time-step restriction (dt ~ h_min / |u| for advection).
    double minGLLSpacing() const;

    // Evaluate the two bases at arbitrary points (used for output/plotting).
    static void evalModifiedBasis(int P, const Vector &z, Matrix &B, Matrix &dB);
    static Matrix evalOrthonormalLegendre(int P, const Vector &z);
};
