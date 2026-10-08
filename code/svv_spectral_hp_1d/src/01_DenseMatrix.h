// =============================================================================
// 01_DenseMatrix.h  --  a deliberately small dense linear-algebra toolkit
//
// Why write our own instead of using Eigen/LAPACK?
//   * Zero dependencies: the project builds with only a C++17 compiler.
//   * Everything that happens to a matrix is visible in ~200 lines of code.
//   * The problems here are tiny (< 500 unknowns), so O(n^3) dense LU is free.
//
// Nektar++ uses `NekMatrix` / `Array<OneD, NekDouble>` for the same job. The
// concepts are identical: row-major storage, an LU factorisation that is
// computed once and reused for every solve.
// =============================================================================
#pragma once

#include <vector>
#include <string>

using Vector = std::vector<double>;

class Matrix
{
public:
    Matrix() = default;
    Matrix(int rows, int cols, double value = 0.0);

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    // Element access: A(i,j). Row-major storage, no bounds checking (speed).
    double &operator()(int i, int j) { return a_[i * cols_ + j]; }
    double operator()(int i, int j) const { return a_[i * cols_ + j]; }

    static Matrix identity(int n);
    static Matrix diagonal(const Vector &d);

    Matrix transpose() const;
    Matrix operator*(const Matrix &B) const;
    Vector operator*(const Vector &x) const;
    Matrix operator+(const Matrix &B) const;
    Matrix operator-(const Matrix &B) const;
    Matrix scaled(double s) const;
    Matrix inverse() const;   // Gauss-Jordan via LU; only for small matrices

    // Sub-matrix with rows [r0,r1) and columns [c0,c1)
    Matrix block(int r0, int r1, int c0, int c1) const;

    double maxAbs() const;    // max_ij |A_ij| -- handy for "is this zero?" tests
    bool isSymmetric(double tol) const;

private:
    int rows_ = 0, cols_ = 0;
    std::vector<double> a_;
};

// LU factorisation with partial pivoting, PA = LU.  Factor once, solve often.
class LUSolver
{
public:
    explicit LUSolver(const Matrix &A);
    Vector solve(const Vector &b) const;
    int size() const { return n_; }

private:
    int n_ = 0;
    Matrix lu_;              // L (unit lower, below diag) and U (upper) packed
    std::vector<int> piv_;   // row permutation
};

// Small vector helpers -------------------------------------------------------
double dot(const Vector &a, const Vector &b);
double norm2(const Vector &a);                   // Euclidean norm
double maxAbs(const Vector &a);
void axpy(double alpha, const Vector &x, Vector &y);   // y += alpha * x
Vector scaled(const Vector &x, double alpha);
Vector hadamard(const Vector &a, const Vector &b);     // element-wise product
