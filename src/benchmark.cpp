#include "nwchem_gpu.h"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv) {
    try {
        if (argc > 4) throw std::invalid_argument("too many arguments");
        const std::string mode = argc > 1 ? argv[1] : "auto";
        const int backend = mode == "auto" ? NWGPU_AUTO : mode == "cpu" ? NWGPU_CPU : mode == "cuda" ? NWGPU_CUDA : -1;
        const int n = argc > 2 ? std::stoi(argv[2]) : 512;
        const int repeats = argc > 3 ? std::stoi(argv[3]) : 3;
        if (backend < 0 || n < 1 || n > 16384 || repeats < 1 || repeats > 10000)
            throw std::invalid_argument("usage: nwgpu_bench [auto|cpu|cuda] [size 1..16384] [repeats 1..10000]");
        const size_t count = static_cast<size_t>(n) * n;
        std::vector<double> a(count), b(count), c(count), reference(count);
        for (size_t i = 0; i < count; ++i) { a[i] = std::sin(double(i)*.01); b[i] = std::cos(double(i)*.02); }
        auto run = [&](int selected, double* output) {
            if (nwgpu_dgemm(selected, 'N', 'N', n, n, n, 1, a.data(), n, b.data(), n, 0, output, n))
                throw std::runtime_error(nwgpu_last_error());
        };
        run(NWGPU_CPU, reference.data());
        run(backend, c.data()); // Warm up allocation and CUDA context outside timing.
        const auto before = nwgpu_get_stats();
        const auto start = std::chrono::steady_clock::now();
        for (int r = 0; r < repeats; ++r) run(backend, c.data());
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()/repeats;
        const auto after = nwgpu_get_stats();
        double error = 0, scale = 0;
        for (size_t i = 0; i < count; ++i) { error = std::max(error, std::abs(c[i]-reference[i])); scale = std::max(scale, std::abs(reference[i])); }
        if (error > 1e-10 * std::max(1., scale)) throw std::runtime_error("CPU/CUDA correctness comparison failed");
        std::cout << std::setprecision(12) << "{\"requested_backend\":\"" << mode << "\",\"size\":" << n
                  << ",\"seconds_per_call_including_transfers\":" << seconds
                  << ",\"gflops\":" << 2.*n*n*n/seconds/1.e9 << ",\"max_abs_error\":" << error
                  << ",\"cpu_calls\":" << after.cpu_gemm-before.cpu_gemm
                  << ",\"cuda_calls\":" << after.cuda_gemm-before.cuda_gemm << "}\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
