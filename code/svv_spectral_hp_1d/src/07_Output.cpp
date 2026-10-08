#include "07_Output.h"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

Samples sampleField(const CGAssembly &A, const Vector &uhat, int nsPerElem)
{
    return sampleField(A, uhat.data(), nsPerElem);
}

Samples sampleField(const CGAssembly &A, const double *uhat, int nsPerElem)
{
    Vector xi(nsPerElem);
    for (int i = 0; i < nsPerElem; ++i)
        xi[i] = -1.0 + 2.0 * i / (nsPerElem - 1);
    Matrix Bs, dBs;
    Basis1D::evalModifiedBasis(A.P, xi, Bs, dBs);

    Samples s;
    for (int e = 0; e < A.nel; ++e)
    {
        const Vector ue = Bs * A.gather(uhat, e);
        for (int i = 0; i < nsPerElem; ++i)
        {
            s.x.push_back(A.mesh.x(e, xi[i]));
            s.u.push_back(ue[i]);
        }
    }
    return s;
}

void writeSpectrumCSV(const std::string &path, const CGAssembly &A, const Vector &uhat)
{
    std::ofstream f(path);
    if (!f)
        throw std::runtime_error("cannot open " + path);
    f << "elem,k,coef\n";
    f << std::setprecision(10);
    for (int e = 0; e < A.nel; ++e)
    {
        const Vector ut = A.basis.toOrthonormal(A.gather(uhat, e));
        for (int k = 0; k < A.nModes; ++k)
            f << e << "," << k << "," << ut[k] << "\n";
    }
}

void writeCSV(const std::string &path, const std::vector<std::pair<std::string, const Vector *>> &cols)
{
    std::ofstream f(path);
    if (!f)
        throw std::runtime_error("cannot open " + path);
    f << std::setprecision(12);
    for (size_t c = 0; c < cols.size(); ++c)
        f << cols[c].first << (c + 1 < cols.size() ? "," : "\n");
    const size_t n = cols.front().second->size();
    for (size_t i = 0; i < n; ++i)
        for (size_t c = 0; c < cols.size(); ++c)
            f << (*cols[c].second)[i] << (c + 1 < cols.size() ? "," : "\n");
}

double totalVariation(const Vector &u)
{
    double tv = 0.0;
    for (size_t i = 1; i < u.size(); ++i)
        tv += std::fabs(u[i] - u[i - 1]);
    return tv;
}

// ---------------------------------------------------------------------------
Args::Args(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i)
        argv_.emplace_back(argv[i]);
}

bool Args::has(const std::string &key) const
{
    for (const auto &a : argv_)
        if (a == key)
            return true;
    return false;
}

std::string Args::getString(const std::string &key, const std::string &def) const
{
    for (size_t i = 0; i + 1 < argv_.size(); ++i)
        if (argv_[i] == key)
            return argv_[i + 1];
    return def;
}

double Args::getDouble(const std::string &key, double def) const
{
    const std::string s = getString(key, "");
    return s.empty() ? def : std::stod(s);
}

int Args::getInt(const std::string &key, int def) const
{
    const std::string s = getString(key, "");
    return s.empty() ? def : std::stoi(s);
}

bool Args::getFlag(const std::string &key) const
{
    for (size_t i = 0; i < argv_.size(); ++i)
        if (argv_[i] == key)
        {
            if (i + 1 < argv_.size())
            {
                const std::string &v = argv_[i + 1];
                if (v == "0" || v == "false" || v == "off")
                    return false;
            }
            return true;
        }
    return false;
}
