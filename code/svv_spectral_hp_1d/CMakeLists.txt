cmake_minimum_required(VERSION 3.16)
project(svv_spectral_hp_1d CXX)

# ---------------------------------------------------------------------------
# Build configuration
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ---------------------------------------------------------------------------
# polylib: Karniadakis & Sherwin's polynomial library (Jacobi polynomials,
# Gauss-Lobatto-Legendre points/weights, derivative matrices). It is the same
# library that lives inside Nektar++ (library/LibUtilities/Polylib). We compile
# it as C++ because polylib.h wraps everything in `namespace polylib`.
# ---------------------------------------------------------------------------
add_library(polylib STATIC extern/polylib/polylib.cpp)
target_include_directories(polylib PUBLIC extern/polylib)
target_compile_options(polylib PRIVATE -w)   # legacy C code, silence warnings

# ---------------------------------------------------------------------------
# svvhp: the spectral/hp + SVV library written for this project.
# Files are numbered in the order you should read them (see README.md).
# ---------------------------------------------------------------------------
add_library(svvhp STATIC
  src/01_DenseMatrix.cpp
  src/02_Basis1D.cpp
  src/03_SVV.cpp
  src/05_CGAssembly.cpp
  src/07_Output.cpp
  src/08_Burgers.cpp
  src/09_CompressibleNS.cpp
)
target_include_directories(svvhp PUBLIC src)
target_link_libraries(svvhp PUBLIC polylib)
target_compile_options(svvhp PRIVATE -Wall -Wextra -O2)

# ---------------------------------------------------------------------------
# Applications
# ---------------------------------------------------------------------------
add_executable(burgers_svv apps/burgers_svv.cpp)
target_link_libraries(burgers_svv PRIVATE svvhp)

add_executable(cns1d_svv apps/cns1d_svv.cpp)
target_link_libraries(cns1d_svv PRIVATE svvhp)

# ---------------------------------------------------------------------------
# Unit tests (run with `ctest --test-dir build --output-on-failure`)
# ---------------------------------------------------------------------------
enable_testing()
add_executable(test_svv_operator tests/test_svv_operator.cpp)
target_link_libraries(test_svv_operator PRIVATE svvhp)
add_test(NAME svv_operator COMMAND test_svv_operator)
