// =============================================================================
// 07_Output.h  --  sampling the modal solution for plots, modal spectra, CSV
// =============================================================================
#pragma once

#include "01_DenseMatrix.h"
#include "05_CGAssembly.h"

#include <string>
#include <utility>
#include <vector>

// Evaluate u_h on nsPerElem uniformly spaced points per element.
struct Samples
{
    Vector x, u;
};
Samples sampleField(const CGAssembly &A, const Vector &uhat, int nsPerElem = 60);
Samples sampleField(const CGAssembly &A, const double *uhat, int nsPerElem = 60);

// Orthonormal (Legendre) coefficients of every element: rows (elem, k, coef).
// This is the "spectrum" the SVV kernel acts on.
void writeSpectrumCSV(const std::string &path, const CGAssembly &A, const Vector &uhat);

// Generic CSV writer: columns are (name, data) pairs of equal length.
void writeCSV(const std::string &path, const std::vector<std::pair<std::string, const Vector *>> &cols);

// Total variation sum |u_{i+1} - u_i| of a sampled signal (a wiggle detector).
double totalVariation(const Vector &u);

// Command-line helper: --key value  (returns default when absent) -------------
class Args
{
public:
    Args(int argc, char **argv);
    bool has(const std::string &key) const;
    double getDouble(const std::string &key, double def) const;
    int getInt(const std::string &key, int def) const;
    std::string getString(const std::string &key, const std::string &def) const;
    bool getFlag(const std::string &key) const;   // present with value 1/true, or bare

private:
    std::vector<std::string> argv_;
};
