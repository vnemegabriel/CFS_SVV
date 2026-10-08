#include "01_DenseMatrix.h"

#include <cmath>
#include <stdexcept>
#include <algorithm>

Matrix::Matrix(int rows, int cols, double value)
    : rows_(rows), cols_(cols), a_(static_cast<size_t>(rows) * cols, value)
{
}

Matrix Matrix::identity(int n)
{
    Matrix I(n, n);
    for (int i = 0; i < n; ++i)
        I(i, i) = 1.0;
    return I;
}

Matrix Matrix::diagonal(const Vector &d)
{
    const int n = static_cast<int>(d.size());
    Matrix D(n, n);
    for (int i = 0; i < n; ++i)
        D(i, i) = d[i];
    return D;
}

Matrix Matrix::transpose() const
{
    Matrix T(cols_, rows_);
    for (int i = 0; i < rows_; ++i)
        for (int j = 0; j < cols_; ++j)
            T(j, i) = (*this)(i, j);
    return T;
}

Matrix Matrix::operator*(const Matrix &B) const
{
    if (cols_ != B.rows_)
        throw std::runtime_error("Matrix*Matrix: dimension mismatch");
    Matrix C(rows_, B.cols_);
    for (int i = 0; i < rows_; ++i)
        for (int k = 0; k < cols_; ++k)
        {
            const double aik = (*this)(i, k);
            if (aik == 0.0)
                continue;
            for (int j = 0; j < B.cols_; ++j)
                C(i, j) += aik * B(k, j);
        }
    return C;
}

Vector Matrix::operator*(const Vector &x) const
{
    if (cols_ != static_cast<int>(x.size()))
        throw std::runtime_error("Matrix*Vector: dimension mismatch");
    Vector y(rows_, 0.0);
    for (int i = 0; i < rows_; ++i)
    {
        double s = 0.0;
        for (int j = 0; j < cols_; ++j)
            s += (*this)(i, j) * x[j];
        y[i] = s;
    }
    return y;
}

Matrix Matrix::operator+(const Matrix &B) const
{
    Matrix C(*this);
    for (size_t i = 0; i < a_.size(); ++i)
        C.a_[i] += B.a_[i];
    return C;
}

Matrix Matrix::operator-(const Matrix &B) const
{
    Matrix C(*this);
    for (size_t i = 0; i < a_.size(); ++i)
        C.a_[i] -= B.a_[i];
    return C;
}

Matrix Matrix::scaled(double s) const
{
    Matrix C(*this);
    for (double &v : C.a_)
        v *= s;
    return C;
}

Matrix Matrix::inverse() const
{
    if (rows_ != cols_)
        throw std::runtime_error("inverse: matrix not square");
    LUSolver lu(*this);
    Matrix Inv(rows_, rows_);
    Vector e(rows_, 0.0);
    for (int j = 0; j < rows_; ++j)
    {
        std::fill(e.begin(), e.end(), 0.0);
        e[j]     = 1.0;
        Vector c = lu.solve(e);
        for (int i = 0; i < rows_; ++i)
            Inv(i, j) = c[i];
    }
    return Inv;
}

Matrix Matrix::block(int r0, int r1, int c0, int c1) const
{
    Matrix B(r1 - r0, c1 - c0);
    for (int i = r0; i < r1; ++i)
        for (int j = c0; j < c1; ++j)
            B(i - r0, j - c0) = (*this)(i, j);
    return B;
}

double Matrix::maxAbs() const
{
    double m = 0.0;
    for (double v : a_)
        m = std::max(m, std::fabs(v));
    return m;
}

bool Matrix::isSymmetric(double tol) const
{
    if (rows_ != cols_)
        return false;
    for (int i = 0; i < rows_; ++i)
        for (int j = i + 1; j < cols_; ++j)
            if (std::fabs((*this)(i, j) - (*this)(j, i)) > tol)
                return false;
    return true;
}

// ---------------------------------------------------------------------------
// LU with partial pivoting (Doolittle). Textbook algorithm, kept readable.
// ---------------------------------------------------------------------------
LUSolver::LUSolver(const Matrix &A) : n_(A.rows()), lu_(A), piv_(A.rows())
{
    if (A.rows() != A.cols())
        throw std::runtime_error("LUSolver: matrix not square");

    for (int i = 0; i < n_; ++i)
        piv_[i] = i;

    for (int k = 0; k < n_; ++k)
    {
        // pivot search in column k
        int p       = k;
        double pmax = std::fabs(lu_(k, k));
        for (int i = k + 1; i < n_; ++i)
            if (std::fabs(lu_(i, k)) > pmax)
            {
                pmax = std::fabs(lu_(i, k));
                p    = i;
            }
        if (pmax == 0.0)
            throw std::runtime_error("LUSolver: singular matrix");
        if (p != k)
        {
            for (int j = 0; j < n_; ++j)
                std::swap(lu_(k, j), lu_(p, j));
            std::swap(piv_[k], piv_[p]);
        }
        // elimination
        for (int i = k + 1; i < n_; ++i)
        {
            lu_(i, k) /= lu_(k, k);
            const double lik = lu_(i, k);
            if (lik == 0.0)
                continue;
            for (int j = k + 1; j < n_; ++j)
                lu_(i, j) -= lik * lu_(k, j);
        }
    }
}

Vector LUSolver::solve(const Vector &b) const
{
    Vector x(n_);
    // forward substitution with L (unit diagonal), applying the permutation
    for (int i = 0; i < n_; ++i)
    {
        double s = b[piv_[i]];
        for (int j = 0; j < i; ++j)
            s -= lu_(i, j) * x[j];
        x[i] = s;
    }
    // back substitution with U
    for (int i = n_ - 1; i >= 0; --i)
    {
        double s = x[i];
        for (int j = i + 1; j < n_; ++j)
            s -= lu_(i, j) * x[j];
        x[i] = s / lu_(i, i);
    }
    return x;
}

// ---------------------------------------------------------------------------
double dot(const Vector &a, const Vector &b)
{
    double s = 0.0;
    for (size_t i = 0; i < a.size(); ++i)
        s += a[i] * b[i];
    return s;
}

double norm2(const Vector &a)
{
    return std::sqrt(dot(a, a));
}

double maxAbs(const Vector &a)
{
    double m = 0.0;
    for (double v : a)
        m = std::max(m, std::fabs(v));
    return m;
}

void axpy(double alpha, const Vector &x, Vector &y)
{
    for (size_t i = 0; i < x.size(); ++i)
        y[i] += alpha * x[i];
}

Vector scaled(const Vector &x, double alpha)
{
    Vector y(x);
    for (double &v : y)
        v *= alpha;
    return y;
}

Vector hadamard(const Vector &a, const Vector &b)
{
    Vector c(a.size());
    for (size_t i = 0; i < a.size(); ++i)
        c[i] = a[i] * b[i];
    return c;
}
