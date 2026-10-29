#include "internal.hpp"
#include "dirac_point.hpp"
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <algorithm>
#include <cstdlib>
#include <climits>
#include <stdexcept>
#include <vector>

namespace nwgpu {
namespace {
void check(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(cudaGetErrorString(status));
}
void check(cublasStatus_t status) {
    if (status != CUBLAS_STATUS_SUCCESS)
        throw std::runtime_error("cuBLAS error " + std::to_string(static_cast<int>(status)));
}
struct Buffer {
    double* ptr = nullptr;
    size_t capacity = 0;
    ~Buffer() { if (ptr) cudaFree(ptr); }
    void reserve(size_t count) {
        if (count <= capacity) return;
        double* next = nullptr;
        check(cudaMalloc(reinterpret_cast<void**>(&next), count * sizeof(double)));
        if (ptr) cudaFree(ptr);
        ptr = next;
        capacity = count;
    }
    Buffer() = default;
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
};
int configured_device() {
    const char* raw = std::getenv("NWGPU_DEVICE");
    if (!raw) return 0;
    char* end = nullptr;
    const long value = std::strtol(raw, &end, 10);
    if (end == raw || *end || value < 0 || value > INT_MAX)
        throw std::invalid_argument("NWGPU_DEVICE must be a nonnegative device ordinal");
    return static_cast<int>(value);
}
struct Context {
    int device;
    cublasHandle_t handle = nullptr;
    Buffer a, b, c, rho, w, a1, a2, a3, f, terms;
    Context() : device(configured_device()) {
        check(cudaSetDevice(device));
        check(cublasCreate(&handle));
        try { check(cublasSetMathMode(handle, CUBLAS_DEFAULT_MATH)); }
        catch (...) { cublasDestroy(handle); throw; }
    }
    ~Context() { cudaSetDevice(device); cublasDestroy(handle); }
};
Context& context() {
    thread_local Context ctx;
    check(cudaSetDevice(ctx.device));
    return ctx;
}
bool trans(char c) { return c == 't' || c == 'T' || c == 'c' || c == 'C'; }
void upload(Buffer& buffer, const double* host, size_t count) {
    buffer.reserve(count);
    check(cudaMemcpy(buffer.ptr, host, count * sizeof(double), cudaMemcpyHostToDevice));
}
std::vector<double> download(const Buffer& buffer, size_t count) {
    std::vector<double> host(count);
    check(cudaMemcpy(host.data(), buffer.ptr, count * sizeof(double), cudaMemcpyDeviceToHost));
    return host;
}
__global__ void dirac_kernel(Dirac d, double* sums) {
    __shared__ double values[256];
    const int lane = threadIdx.x;
    const size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + lane;
    values[lane] = i < static_cast<size_t>(d.nq) ? dirac_point(d, static_cast<int>(i)) : 0.;
    __syncthreads();
    for (int offset = 128; offset > 0; offset /= 2) {
        if (lane < offset) values[lane] += values[lane + offset];
        __syncthreads();
    }
    if (lane == 0) sums[blockIdx.x] = values[0];
}
}
bool cuda_available() {
    try {
        int count = 0;
        if (cudaGetDeviceCount(&count) != cudaSuccess) { cudaGetLastError(); return false; }
        return configured_device() < count;
    } catch (...) { return false; }
}
void cuda_gemm(const Gemm& g) {
    auto& ctx = context();
    const int ar = trans(g.ta) ? g.k : g.m, ac = trans(g.ta) ? g.m : g.k;
    const int br = trans(g.tb) ? g.n : g.k, bc = trans(g.tb) ? g.k : g.n;
    const size_t bytes = sizeof(double);
    ctx.a.reserve(static_cast<size_t>(ar) * ac);
    ctx.b.reserve(static_cast<size_t>(br) * bc);
    ctx.c.reserve(static_cast<size_t>(g.m) * g.n);
    check(cudaMemcpy2D(ctx.a.ptr, ar * bytes, g.a, g.lda * bytes, ar * bytes, ac, cudaMemcpyHostToDevice));
    check(cudaMemcpy2D(ctx.b.ptr, br * bytes, g.b, g.ldb * bytes, br * bytes, bc, cudaMemcpyHostToDevice));
    if (g.beta != 0)
        check(cudaMemcpy2D(ctx.c.ptr, g.m * bytes, g.c, g.ldc * bytes, g.m * bytes, g.n, cudaMemcpyHostToDevice));
    check(cublasDgemm(ctx.handle, trans(g.ta) ? CUBLAS_OP_T : CUBLAS_OP_N,
                     trans(g.tb) ? CUBLAS_OP_T : CUBLAS_OP_N, g.m, g.n, g.k,
                     &g.alpha, ctx.a.ptr, ar, ctx.b.ptr, br, &g.beta, ctx.c.ptr, g.m));
    // Stage the entire result before committing. Auto fallback sees unchanged C
    // if allocation, launch, or transfer fails, even for beta != 0.
    auto result = download(ctx.c, static_cast<size_t>(g.m) * g.n);
    for (int j = 0; j < g.n; ++j)
        std::copy_n(result.data() + static_cast<size_t>(j) * g.m, g.m,
                    g.c + static_cast<size_t>(j) * g.ldc);
}
void cuda_dirac(const Dirac& d) {
    auto& ctx = context();
    const size_t q = static_cast<size_t>(d.nq);
    const size_t blocks = (q + 255) / 256;
    upload(ctx.rho, d.rho, q * (d.ipol == 1 ? 1 : 3));
    upload(ctx.w, d.w, q);
    upload(ctx.a1, d.a, q * d.ipol);
    if (d.order >= 2) upload(ctx.a2, d.a2, q * 3);
    if (d.order == 3) upload(ctx.a3, d.a3, q * 4);
    if (d.ldew) upload(ctx.f, d.func, q);
    ctx.terms.reserve(blocks);
    Dirac device = d;
    device.rho = ctx.rho.ptr; device.w = ctx.w.ptr;
    device.a = ctx.a1.ptr; device.a2 = ctx.a2.ptr; device.a3 = ctx.a3.ptr;
    device.func = ctx.f.ptr; device.ex = nullptr;
    dirac_kernel<<<static_cast<unsigned>(blocks), 256>>>(device, ctx.terms.ptr);
    check(cudaGetLastError());
    auto sums = download(ctx.terms, blocks);
    auto a1 = download(ctx.a1, q * d.ipol);
    auto a2 = d.order >= 2 ? download(ctx.a2, q * 3) : std::vector<double>{};
    auto a3 = d.order == 3 ? download(ctx.a3, q * 4) : std::vector<double>{};
    auto f = d.ldew ? download(ctx.f, q) : std::vector<double>{};
    double term = 0;
    for (double value : sums) term += value;
    std::copy(a1.begin(), a1.end(), d.a);
    if (d.order >= 2) std::copy(a2.begin(), a2.end(), d.a2);
    if (d.order == 3) std::copy(a3.begin(), a3.end(), d.a3);
    if (d.ldew) std::copy(f.begin(), f.end(), d.func);
    *d.ex += dirac_energy_factor(d) * term;
}
}
