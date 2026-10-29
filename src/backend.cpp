#include "internal.hpp"
#include "dirac_point.hpp"
#include <algorithm>
#include <atomic>
#include <climits>
#include <cstdlib>
#include <cstdio>
#include <limits>
#include <stdexcept>

#ifdef NWGPU_HAVE_BLAS
// LP64 CBLAS avoids binding our CPU path to NWChem's Fortran dgemm_ symbol,
// which may use a different INTEGER width. These are the standard CBLAS enums.
extern "C" void cblas_dgemm(int, int, int, int, int, int, double,
                            const double*, int, const double*, int, double, double*, int);
#endif
namespace {
thread_local std::string error;
std::atomic<uint64_t> cpu_gemm{0}, gpu_gemm{0}, cpu_dirac{0}, gpu_dirac{0}, fallbacks{0};
struct Report {
    ~Report() {
        const char* enabled = std::getenv("NWGPU_REPORT_STATS");
        if (!enabled || std::string(enabled) != "1") return;
        const char* rank = std::getenv("OMPI_COMM_WORLD_RANK");
        if (!rank) rank = std::getenv("PMI_RANK");
        std::fprintf(stderr, "nwchem_gpu rank=%s cpu_gemm=%llu cuda_gemm=%llu cpu_dirac=%llu cuda_dirac=%llu fallbacks=%llu\n",
            rank ? rank : "unknown", static_cast<unsigned long long>(cpu_gemm.load()),
            static_cast<unsigned long long>(gpu_gemm.load()), static_cast<unsigned long long>(cpu_dirac.load()),
            static_cast<unsigned long long>(gpu_dirac.load()), static_cast<unsigned long long>(fallbacks.load()));
    }
} report;
bool transpose(char c) { return c == 't' || c == 'T' || c == 'c' || c == 'C'; }
bool valid_op(char c) { return c == 'n' || c == 'N' || transpose(c); }
void require(bool ok, const char* message) { if (!ok) throw std::invalid_argument(message); }
double threshold(const char* name, double fallback) {
    const char* raw = std::getenv(name);
    if (!raw) return fallback;
    char* end = nullptr;
    const double value = std::strtod(raw, &end);
    require(end != raw && *end == '\0' && std::isfinite(value) && value >= 0,
            "Invalid nonnegative dispatch threshold");
    return value;
}
void cpu_multiply(const nwgpu::Gemm& g) {
    // Quick paths must not access A/B when alpha=0 or k=0, or C when beta=0.
    if (g.alpha == 0 || g.k == 0) {
        if (g.beta == 1) return;
        for (int j = 0; j < g.n; ++j)
            for (int i = 0; i < g.m; ++i) {
                double& c = g.c[static_cast<size_t>(j) * g.ldc + i];
                c = g.beta == 0 ? 0 : g.beta * c;
            }
        return;
    }
#ifdef NWGPU_HAVE_BLAS
    cblas_dgemm(102, transpose(g.ta) ? 112 : 111, transpose(g.tb) ? 112 : 111,
                g.m, g.n, g.k, g.alpha, g.a, g.lda, g.b, g.ldb, g.beta, g.c, g.ldc);
#else
    for (int j = 0; j < g.n; ++j)
        for (int i = 0; i < g.m; ++i) {
            double sum = 0;
            for (int p = 0; p < g.k; ++p) {
                const auto ai = transpose(g.ta) ? static_cast<size_t>(i)*g.lda+p : static_cast<size_t>(p)*g.lda+i;
                const auto bi = transpose(g.tb) ? static_cast<size_t>(p)*g.ldb+j : static_cast<size_t>(j)*g.ldb+p;
                sum += g.a[ai] * g.b[bi];
            }
            double& c = g.c[static_cast<size_t>(j)*g.ldc+i];
            c = g.alpha * sum + (g.beta == 0 ? 0 : g.beta * c);
        }
#endif
}
bool use_cuda(int backend, bool large) {
    require(backend >= NWGPU_AUTO && backend <= NWGPU_CUDA, "Invalid backend");
    if (backend == NWGPU_CPU) return false;
    if (backend == NWGPU_CUDA) {
        require(nwgpu_cuda_available(), "CUDA requested but no usable CUDA backend/device is available");
        return true;
    }
    return large && nwgpu_cuda_available();
}
}
namespace nwgpu {
void fail(const std::string& message) { error = message; }
int env_backend() {
    const char* value = std::getenv("NWGPU_BACKEND");
    if (!value || std::string(value) == "auto") return NWGPU_AUTO;
    if (std::string(value) == "cpu") return NWGPU_CPU;
    if (std::string(value) == "cuda") return NWGPU_CUDA;
    throw std::invalid_argument("NWGPU_BACKEND must be auto, cpu, or cuda");
}
}
extern "C" const char* nwgpu_last_error() { return error.c_str(); }
extern "C" int nwgpu_cuda_available() {
#ifdef NWGPU_HAVE_CUDA
    return nwgpu::cuda_available();
#else
    return 0;
#endif
}
extern "C" nwgpu_stats nwgpu_get_stats() {
    return {cpu_gemm.load(), gpu_gemm.load(), cpu_dirac.load(), gpu_dirac.load(), fallbacks.load()};
}
extern "C" int nwgpu_dgemm(int backend, char ta, char tb, int64_t m, int64_t n, int64_t k,
                            double alpha, const double* a, int64_t lda, const double* b,
                            int64_t ldb, double beta, double* c, int64_t ldc) {
    error.clear();
    try {
        require(valid_op(ta) && valid_op(tb), "DGEMM transpose must be N/T/C");
        require(m >= 0 && n >= 0 && k >= 0 && m <= INT_MAX && n <= INT_MAX && k <= INT_MAX,
                "DGEMM dimensions must be in [0, INT_MAX]");
        require(lda >= std::max<int64_t>(1, transpose(ta) ? k : m) && lda <= INT_MAX, "Invalid DGEMM lda");
        require(ldb >= std::max<int64_t>(1, transpose(tb) ? n : k) && ldb <= INT_MAX, "Invalid DGEMM ldb");
        require(ldc >= std::max<int64_t>(1, m) && ldc <= INT_MAX, "Invalid DGEMM ldc");
        require(backend >= 0 && backend <= 2, "Invalid backend");
        if (!m || !n) return 0;
        if ((alpha == 0 || k == 0) && beta == 1) return 0;
        const auto max_elements = static_cast<uint64_t>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(double);
        auto fits = [max_elements](int64_t ld, int64_t columns) {
            return static_cast<uint64_t>(ld) * static_cast<uint64_t>(columns) <= max_elements;
        };
        require(fits(ldc, n), "DGEMM C extent exceeds addressable memory");
        if (alpha != 0 && k != 0)
            require(fits(lda, transpose(ta) ? m : k) && fits(ldb, transpose(tb) ? k : n),
                    "DGEMM A/B extent exceeds addressable memory");
        require(c, "DGEMM C is null");
        require((alpha == 0 || k == 0) || (a && b), "DGEMM A/B is null");
        nwgpu::Gemm g{ta, tb, int(m), int(n), int(k), int(lda), int(ldb), int(ldc), alpha, beta, a, b, c};
        const bool nontrivial = alpha != 0 && k != 0;
        if (use_cuda(backend, nontrivial && 2. * m * n * k >= threshold("NWGPU_MIN_GEMM_FLOPS", 1.e8))) {
#ifdef NWGPU_HAVE_CUDA
            if (nontrivial) {
                try { nwgpu::cuda_gemm(g); ++gpu_gemm; return 0; }
                catch (...) { if (backend == NWGPU_CUDA) throw; ++fallbacks; }
            }
#endif
        }
        cpu_multiply(g); ++cpu_gemm;
        return 0;
    } catch (const std::exception& e) { nwgpu::fail(e.what()); return 1; }
}
extern "C" int nwgpu_dirac(int backend, int order, int64_t nq, int ipol, double tol,
                            double fac, const double* rho, const double* weights,
                            double* amat, double* amat2, double* amat3, double* ex,
                            int energy_density, double* func) {
    error.clear();
    try {
        require(order >= 1 && order <= 3 && (ipol == 1 || ipol == 2), "Invalid Dirac order/spin");
        require(nq >= 0 && nq <= INT_MAX, "Dirac nq must be in [0, INT_MAX]");
        require(std::isfinite(tol) && tol >= 0, "Dirac density tolerance must be finite and nonnegative");
        require(backend >= 0 && backend <= 2, "Invalid backend");
        if (!nq) return 0;
        require(rho && weights && amat && ex, "Null Dirac input/output");
        require(order < 2 || amat2, "Dirac second derivatives require amat2");
        require(order < 3 || amat3, "Dirac third derivatives require amat3");
        require(!energy_density || func, "Dirac energy density requires func");
        nwgpu::Dirac d{order, int(nq), ipol, tol, fac, rho, weights, amat, amat2, amat3, ex, energy_density != 0, func};
        if (use_cuda(backend, nq >= threshold("NWGPU_MIN_GRID_POINTS", 32768))) {
#ifdef NWGPU_HAVE_CUDA
            try { nwgpu::cuda_dirac(d); ++gpu_dirac; return 0; }
            catch (...) { if (backend == NWGPU_CUDA) throw; ++fallbacks; }
#endif
        }
        double term = 0;
        for (int i = 0; i < d.nq; ++i) term += nwgpu::dirac_point(d, i);
        *ex += nwgpu::dirac_energy_factor(d) * term;
        ++cpu_dirac;
        return 0;
    } catch (const std::exception& e) { nwgpu::fail(e.what()); return 1; }
}
