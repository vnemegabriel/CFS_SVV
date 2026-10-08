// =============================================================================
// 04_Mesh1D.h  --  a 1D mesh is just a sorted list of vertices.
//
// Element e occupies [xv[e], xv[e+1]].  The affine map from the reference
// interval is
//        x(xi) = xv[e] + (1+xi)/2 * h_e ,   dx/dxi = h_e/2 =: J_e   (Jacobian)
// so that   d/dx = (1/J_e) d/dxi   and   dx = J_e dxi.
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"

#include <stdexcept>

struct Mesh1D
{
    Vector xv;   // vertex coordinates, size nel+1, increasing

    int nel() const { return static_cast<int>(xv.size()) - 1; }
    double h(int e) const { return xv[e + 1] - xv[e]; }
    double J(int e) const { return 0.5 * h(e); }
    double x(int e, double xi) const { return xv[e] + 0.5 * (1.0 + xi) * h(e); }
    double xmin() const { return xv.front(); }
    double xmax() const { return xv.back(); }

    static Mesh1D uniform(double a, double b, int nel)
    {
        if (nel < 1 || b <= a)
            throw std::runtime_error("Mesh1D::uniform: bad arguments");
        Mesh1D m;
        m.xv.resize(nel + 1);
        for (int i = 0; i <= nel; ++i)
            m.xv[i] = a + (b - a) * i / nel;
        return m;
    }
};
