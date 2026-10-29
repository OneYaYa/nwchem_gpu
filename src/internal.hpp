#pragma once
#include "nwchem_gpu.h"
#include <cmath>
#include <cstddef>
#include <string>

namespace nwgpu {
struct Gemm {
    char ta, tb;
    int m, n, k, lda, ldb, ldc;
    double alpha, beta;
    const double *a, *b;
    double* c;
};
struct Dirac {
    int order, nq, ipol;
    double tol, fac;
    const double *rho, *w;
    double *a, *a2, *a3, *ex;
    bool ldew;
    double* func;
};
int env_backend();
void fail(const std::string& message);
#ifdef NWGPU_HAVE_CUDA
bool cuda_available();
void cuda_gemm(const Gemm& g);
void cuda_dirac(const Dirac& d);
#endif
}
