#include "internal.hpp"
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>

#if NWGPU_FORTRAN_INTEGER_BYTES == 8
using FInt = int64_t;
#else
using FInt = int32_t;
#endif
namespace {
[[noreturn]] void fatal(const char* routine) {
    std::fprintf(stderr, "nwchem_gpu %s: %s\n", routine, nwgpu_last_error());
    std::abort();
}
int backend() {
    try { return nwgpu::env_backend(); }
    catch (const std::exception& e) { nwgpu::fail(e.what()); fatal("configuration"); }
}
void dirac(int order, const double* tol, const double* fac, const double* rho,
           double* a, double* a2, double* a3, const FInt* nq, const FInt* ipol,
           double* ex, const double* w, const FInt* ldew, double* func) {
    if (*ipol != 1 && *ipol != 2) { nwgpu::fail("Invalid ipol"); fatal("xc_dirac"); }
    if (nwgpu_dirac(backend(), order, *nq, static_cast<int>(*ipol), *tol, *fac,
                    rho, w, a, a2, a3, ex, *ldew != 0, func)) fatal("xc_dirac");
}
}
// nwgmm has the same length as dgemm, preserving fixed-form Fortran columns.
// These arguments are NWChem INTEGERs, independent of BLAS_SIZE/USE_64TO32.
extern "C" void nwgmm_(const char* ta, const char* tb, const FInt* m,
                         const FInt* n, const FInt* k, const double* alpha,
                         const double* a, const FInt* lda, const double* b,
                         const FInt* ldb, const double* beta, double* c, const FInt* ldc) {
    if (nwgpu_dgemm(backend(), *ta, *tb, *m, *n, *k, *alpha, a, *lda, b, *ldb,
                    *beta, c, *ldc)) fatal("nwgmm");
}
extern "C" void xc_dirac_(const double* tol, const double* fac, const FInt*, const FInt*,
                           const double* rho, double* a, const FInt* nq, const FInt* ipol,
                           double* ex, const double* w, const FInt* ldew, double* func) {
    dirac(1, tol, fac, rho, a, nullptr, nullptr, nq, ipol, ex, w, ldew, func);
}
extern "C" void xc_dirac_d2_(const double* tol, const double* fac, const FInt*, const FInt*,
                              const double* rho, double* a, double* a2, const FInt* nq,
                              const FInt* ipol, double* ex, const double* w,
                              const FInt* ldew, double* func) {
    dirac(2, tol, fac, rho, a, a2, nullptr, nq, ipol, ex, w, ldew, func);
}
extern "C" void xc_dirac_d3_(const double* tol, const double* fac, const FInt*, const FInt*,
                              const double* rho, double* a, double* a2, double* a3,
                              const FInt* nq, const FInt* ipol, double* ex,
                              const double* w, const FInt* ldew, double* func) {
    dirac(3, tol, fac, rho, a, a2, a3, nq, ipol, ex, w, ldew, func);
}
