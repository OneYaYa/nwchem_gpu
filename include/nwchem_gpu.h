#pragma once
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif
/* Host arrays, double precision, column major. Calls are synchronous.
 * Backend: 0=automatic, 1=CPU, 2=CUDA required. No silent CPU fallback in 2.
 * Return 0 on success; error text is local to the calling thread. */
enum nwgpu_backend { NWGPU_AUTO = 0, NWGPU_CPU = 1, NWGPU_CUDA = 2 };
struct nwgpu_stats { uint64_t cpu_gemm, cuda_gemm, cpu_dirac, cuda_dirac, fallbacks; };
int nwgpu_cuda_available(void);
const char* nwgpu_last_error(void);
struct nwgpu_stats nwgpu_get_stats(void);
int nwgpu_dgemm(int backend, char ta, char tb, int64_t m, int64_t n, int64_t k,
                double alpha, const double* a, int64_t lda, const double* b,
                int64_t ldb, double beta, double* c, int64_t ldc);
/* Exact xc_dirac.F accumulation conventions, including restricted rescaling.
 * ipol=1: rho[nq]; ipol=2: rho[nq,3] = total, alpha, beta.
 * amat[nq,ipol], amat2[nq,3], amat3[nq,4]. Unused columns are preserved.
 * order=1,2,3; amat2/amat3 may be null below their respective order.
 * func may be null when energy_density=0. Ex is accumulated, not replaced. */
int nwgpu_dirac(int backend, int order, int64_t nq, int ipol, double tol,
                double fac, const double* rho, const double* weights,
                double* amat, double* amat2, double* amat3, double* ex,
                int energy_density, double* func);
#ifdef __cplusplus
}
#endif
