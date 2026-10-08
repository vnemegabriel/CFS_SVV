#include "05_CGAssembly.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

CGAssembly::CGAssembly(const Mesh1D &mesh_, const Basis1D &basis_, BCType bc_)
    : mesh(mesh_), basis(basis_), bc(bc_), nel(mesh_.nel()), P(basis_.P), nModes(basis_.nModes)
{
    nGlobal = (bc == BCType::Periodic) ? nel * P : nel * P + 1;
}

int CGAssembly::dof(int e, int p) const
{
    const int g = e * P + p;
    return (bc == BCType::Periodic) ? (g % nGlobal) : g;
}

Vector CGAssembly::gather(const Vector &global, int e) const
{
    return gather(global.data(), e);
}

Vector CGAssembly::gather(const double *global, int e) const
{
    Vector local(nModes);
    for (int p = 0; p < nModes; ++p)
        local[p] = global[dof(e, p)];
    return local;
}

void CGAssembly::scatterAdd(const Vector &local, int e, Vector &global) const
{
    scatterAdd(local, e, global.data());
}

void CGAssembly::scatterAdd(const Vector &local, int e, double *global) const
{
    for (int p = 0; p < nModes; ++p)
        global[dof(e, p)] += local[p];
}

Matrix CGAssembly::assemble(const std::function<Matrix(int)> &elementMatrix) const
{
    Matrix A(nGlobal, nGlobal);
    for (int e = 0; e < nel; ++e)
    {
        const Matrix Ae = elementMatrix(e);
        for (int i = 0; i < nModes; ++i)
            for (int j = 0; j < nModes; ++j)
                A(dof(e, i), dof(e, j)) += Ae(i, j);
    }
    return A;
}

Matrix CGAssembly::massMatrix() const
{
    return assemble([this](int e) { return basis.M.scaled(mesh.J(e)); });
}

Matrix CGAssembly::laplacianMatrix() const
{
    return assemble([this](int e) { return basis.K.scaled(1.0 / mesh.J(e)); });
}

Matrix CGAssembly::operatorMatrix(const Matrix &Lref) const
{
    return assemble([this, &Lref](int e) { return Lref.scaled(1.0 / mesh.J(e)); });
}

std::vector<int> CGAssembly::boundaryDofs() const
{
    if (bc == BCType::Periodic)
        return {};
    return {0, nel * P};
}

void CGAssembly::applyDirichletRows(Matrix &A) const
{
    for (int d : boundaryDofs())
    {
        for (int j = 0; j < nGlobal; ++j)
            A(d, j) = 0.0;
        A(d, d) = 1.0;
    }
}

void CGAssembly::zeroDirichlet(Vector &r) const
{
    zeroDirichlet(r.data());
}

void CGAssembly::zeroDirichlet(double *r) const
{
    for (int d : boundaryDofs())
        r[d] = 0.0;
}

Vector CGAssembly::l2Projection(const std::function<double(double)> &f) const
{
    // Solve  M_global u^ = sum_e J_e (psi_i, f)_e
    Vector rhs(nGlobal, 0.0);
    Vector fq(basis.Q);
    for (int e = 0; e < nel; ++e)
    {
        for (int q = 0; q < basis.Q; ++q)
            fq[q] = f(mesh.x(e, basis.zq[q]));
        scatterAdd(scaled(basis.innerProduct(fq), mesh.J(e)), e, rhs);
    }
    LUSolver lu(massMatrix());
    return lu.solve(rhs);
}

double CGAssembly::l2Norm(const Vector &uhat, const std::function<double(double)> *ref) const
{
    double sum = 0.0;
    for (int e = 0; e < nel; ++e)
    {
        const Vector uq = basis.evaluate(gather(uhat, e));
        for (int q = 0; q < basis.Q; ++q)
        {
            double d = uq[q];
            if (ref)
                d -= (*ref)(mesh.x(e, basis.zq[q]));
            sum += basis.wq[q] * mesh.J(e) * d * d;
        }
    }
    return std::sqrt(sum);
}

double CGAssembly::minSpacing() const
{
    double hmin = 1e300;
    for (int e = 0; e < nel; ++e)
        hmin = std::min(hmin, mesh.J(e) * basis.minGLLSpacing());
    return hmin;
}

double CGAssembly::spectralRadius(const Matrix &A, int iterations) const
{
    Matrix Mbc = massMatrix();
    applyDirichletRows(Mbc);
    LUSolver lu(Mbc);

    // power iteration on  v <- M^{-1} A v  (Dirichlet rows of A v zeroed)
    Vector v(nGlobal);
    for (int i = 0; i < nGlobal; ++i)
        v[i] = std::sin(1.0 + 0.7 * i);   // deterministic "random" start
    zeroDirichlet(v);
    double lambda = 0.0;
    for (int it = 0; it < iterations; ++it)
    {
        Vector Av = A * v;
        zeroDirichlet(Av);
        Vector w   = lu.solve(Av);
        const double nw = norm2(w);
        if (nw == 0.0)
            return 0.0;
        lambda = nw / norm2(v);
        v      = scaled(w, 1.0 / nw);
    }
    return lambda;
}
