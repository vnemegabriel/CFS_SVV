// =============================================================================
// 05_CGAssembly.h  --  from local (element) to global (C0-continuous) unknowns
//
// The continuous Galerkin method requires u to be continuous across element
// boundaries. With the modified basis this is trivial: the only modes that are
// non-zero at a vertex are the two vertex modes, so continuity means the right
// vertex coefficient of element e and the left vertex coefficient of element
// e+1 are THE SAME global unknown.
//
// Local-to-global map (the "scatter/gather" of K&S chapter 4):
//        global_dof(e, p) = e*P + p ,  p = 0..P
// Element e:   [ e*P , e*P+1 , ... , e*P+P ]
// Element e+1: [ (e+1)*P = e*P+P , ... ]       <- shared vertex.
// Total unknowns: nel*P + 1 (Dirichlet/open ends) or nel*P (periodic, where
// the last vertex is identified with the first).
//
// Global matrices are assembled by scatter-adding element matrices; global
// vectors of "inner products with test functions" likewise. The global mass
// matrix is banded, but for these problem sizes we keep it dense and LU it once.
//
// BOUNDARY CONDITIONS
//   * Dirichlet: the time derivative of the boundary unknown is zero. We
//     replace the corresponding ROW of the mass matrix by a unit row and set the
//     corresponding RHS entry to 0. This is exact: the unknown's column still
//     multiplies a known zero rate.
//   * Periodic: handled entirely by the dof map (no boundary dofs at all).
//
// NOTE ON INTERFACE FLUX TERMS: when a second-order term is integrated by parts,
// the boundary terms [v G] at element interfaces are NOT added (the global test
// function is continuous; the discrete weak form only has the domain-boundary
// term, which vanishes for Dirichlet rows and is dropped for periodic cases).
// This is the standard C0 Galerkin treatment.
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"
#include "02_Basis1D.h"
#include "04_Mesh1D.h"

#include <functional>
#include <vector>

enum class BCType
{
    Dirichlet,
    Periodic
};

class CGAssembly
{
public:
    CGAssembly(const Mesh1D &mesh, const Basis1D &basis, BCType bc);

    const Mesh1D &mesh;
    const Basis1D &basis;
    BCType bc;
    int nel, P, nModes, nGlobal;

    int dof(int e, int p) const;

    // element <-> global -----------------------------------------------------
    Vector gather(const Vector &global, int e) const;
    Vector gather(const double *global, int e) const;
    void scatterAdd(const Vector &local, int e, Vector &global) const;
    void scatterAdd(const Vector &local, int e, double *global) const;

    // Assemble a global matrix from a per-element matrix (nModes x nModes).
    Matrix assemble(const std::function<Matrix(int)> &elementMatrix) const;

    Matrix massMatrix() const;                      // sum_e  J_e M_ref
    Matrix laplacianMatrix() const;                 // sum_e (1/J_e) K_ref
    Matrix operatorMatrix(const Matrix &Lref) const; // sum_e (1/J_e) Lref  (e.g. SVV)

    // Boundary conditions ------------------------------------------------------
    std::vector<int> boundaryDofs() const;          // {} for periodic
    void applyDirichletRows(Matrix &A) const;
    void zeroDirichlet(Vector &r) const;
    void zeroDirichlet(double *r) const;

    // Global L2 projection of a function f(x) onto the C0 space.
    Vector l2Projection(const std::function<double(double)> &f) const;

    // ||u_h - ref||_L2 over the domain (ref may be null -> ||u_h||).
    double l2Norm(const Vector &uhat, const std::function<double(double)> *ref = nullptr) const;

    // Smallest physical GLL spacing over all elements (for dt estimates).
    double minSpacing() const;

    // Spectral radius of  M_bc^{-1} A_bc  by power iteration. Used to size the
    // explicit time step of a linear operator A (Laplacian, SVV, ...).
    double spectralRadius(const Matrix &A, int iterations = 300) const;
};
