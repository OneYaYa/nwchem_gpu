#!/usr/bin/env bash
set -euo pipefail
project_dir="$(cd "$(dirname "$0")/.." && pwd)"
cd "$project_dir"
build_dir="${1:-build-cuda}"
cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
  -DNWGPU_ENABLE_CUDA=ON -DNWGPU_FORTRAN_INTEGER_BYTES=8
cmake --build "$build_dir" --parallel
ctest --test-dir "$build_dir" --output-on-failure
if command -v compute-sanitizer >/dev/null 2>&1; then
  compute-sanitizer --tool memcheck --error-exitcode 1 "$build_dir/nwgpu_tests" cuda
fi
"$build_dir/nwgpu_bench" cpu 2048 5
"$build_dir/nwgpu_bench" cuda 2048 5
# Set FC to gfortran to compare the original source on the actual GPU too.
if [[ -n "${FC:-}" ]]; then
  NWGPU_MIN_GEMM_FLOPS=0 NWGPU_MIN_GRID_POINTS=0 \
    python3 tests/compare_fortran.py --fc "$FC" \
    --library "$build_dir/libnwchem_gpu.so" --backend cuda --work "$build_dir/oracle"
fi
