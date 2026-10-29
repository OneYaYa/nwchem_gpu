# NWChem C++ / CUDA migration

The goal is GPU acceleration across the original NWChem application. This directory
contains an initial, integrated migration stage, **not a complete C++ rewrite**.
The full original application still needs Fortran, MPI, Global Arrays, BLAS and
LAPACK. No end-to-end speedup has been measured, and the modified full NWChem
executable has not yet been built or validated.

All work is contained here. The original `../src`, inputs, and existing tutorials
are unchanged.

| Location | Contents |
| --- | --- |
| `upstream/` | Separate snapshot of the original tracked application, with GPU call sites connected |
| `src/` | C++17 implementation, Fortran ABI bridge, CUDA kernels and cuBLAS backend |
| `include/nwchem_gpu.h` | Public C/C++ API with explicit backend and error status |
| `migration_manifest.json` | Source revision, hashes, changed files, counts for every source directory |
| `MODULE_STATUS.md` | Complete directory inventory and remaining migration work |
| `tests/` | Numerical contracts, original Fortran comparisons, integration checks |
| `validation/` | Results actually observed on this Mac |

What is implemented:

- Double-precision `DGEMM`: all N/T/C combinations, column-major storage, leading
  dimensions, alpha/beta semantics, and a CPU implementation using LP64 CBLAS.
  CUDA uses `cublasDgemm`, transfers active matrix regions, and reuses device buffers.
- Dirac exchange: restricted and unrestricted spins, energy, potential, second
  and third derivative arrays, density screening, and optional energy density.
  The port preserves the original restricted accumulation/rescaling convention.
- An integration tool that copies every tracked module and routes **13,140 direct
  GEMM calls in 7,446 files across 31 source directories** to `nwgmm_`. All **87**
  source directories and top-level build files are retained. Counts describe source sites, not
  executed calls or runtime coverage. Many are generated coupled-cluster sources.
- `xc_dirac_`, `xc_dirac_d2_`, and `xc_dirac_d3_` resolve to C++/CUDA. The old
  implementations remain in the copied tree as `xc_origx*` for reference.

Only these operations have been ported. Integrals, most exchange-correlation
functionals (including the remaining PBE/B3LYP terms), gradients, FFTs, SCF
control flow, geometry optimization, vibrations, input parsing and distributed
communication remain in their original implementation. Modules with no routed
calls receive no new acceleration from this stage. Existing upstream CUDA code
and build options are retained separately.

## Build and test on this Mac

CUDA cannot run on the available Apple GPU. The local build uses Apple Accelerate
for CPU GEMM and compiles the same C++ exchange implementation used by CUDA.
CMake and a separate Fortran validation compiler were installed in ignored
`build-tools/` and `build-fortran/` directories, without changing the existing
NWChem Conda environment.

```sh
cd /Users/yuepan/Desktop/nwchem/nwchem_gpu
build-tools/bin/cmake -S . -B build-cpu -DCMAKE_BUILD_TYPE=Release
build-tools/bin/cmake --build build-cpu --parallel
build-tools/bin/ctest --test-dir build-cpu --output-on-failure
build-cpu/nwgpu_bench cpu 512 3
```

To repeat the direct original-Fortran comparison with the local compiler:

```sh
/opt/anaconda3/bin/conda run --prefix ./build-fortran \
  /usr/bin/python3 tests/compare_fortran.py \
  --fc ./build-fortran/bin/arm64-apple-darwin20.0.0-gfortran \
  --library ./build-cpu/libnwchem_gpu.dylib
```

The oracle compiles `../src/nwdft/xc/xc_dirac.F` and
`../src/blas/double/dgemm.f` directly, then compares randomized and edge-case
outputs through the actual Fortran-facing symbols. It does not require NumPy.

## Build on an NVIDIA Linux machine

Use an NVIDIA driver and compatible CUDA toolkit with `nvcc`, CMake 3.24+, a
C++17 compiler, and an LP64 CBLAS provider such as OpenBLAS. The full application
also requires the original NWChem build dependencies. Consult the official
[CUDA Linux installation guide](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/)
and [cuBLAS documentation](https://docs.nvidia.com/cuda/cublas/).

```sh
cmake -S . -B build-cuda -DCMAKE_BUILD_TYPE=Release \
  -DNWGPU_ENABLE_CUDA=ON -DNWGPU_FORTRAN_INTEGER_BYTES=8
cmake --build build-cuda --parallel
ctest --test-dir build-cuda --output-on-failure
./build-cuda/nwgpu_bench cpu 2048 5
./build-cuda/nwgpu_bench cuda 2048 5
```

Set `CMAKE_CUDA_ARCHITECTURES` explicitly if building for a GPU other than the
compiler's default target. CUDA is opt-in: requesting it without a toolkit fails
configuration. The CUDA CTest requires a usable GPU and fails if none is present.
`bash tools/validate_cuda.sh` runs these checks and compute-sanitizer when available.
With `FC=gfortran`, it also compares the original Fortran routines on CUDA.
**CUDA compilation, execution, sanitizer results, and performance are unverified
in this workspace.**

## Build the complete application with the new backend

The prepared `upstream/` tree is already present here and points to the local
CPU backend. On another machine, build its backend first and either prepare a
fresh tree or override `NWGPU_LIBRARY` after sourcing the environment file.

```sh
# Only needed to create a fresh tree; the tool refuses to overwrite one.
python3 tools/prepare_nwchem.py --destination upstream-cuda \
  --library "$PWD/build-cuda/libnwchem_gpu.so"

source upstream-cuda/nwchem_gpu_env.sh
# Set your normal NWChem environment, MPI, BLAS and LAPACK options per INSTALL.
# NWCHEM_TOP must continue to point to this copied tree.
# Select all desired NWCHEM_MODULES, including optional modules as before.
cd "$NWCHEM_TOP/src"
make nwchem_config
# If the normal build uses USE_64TO32, run its usual make 64_to_32 step here.
make -j4
```

The snapshot includes tracked working-tree edits, the original build files,
scientific data and QA fixtures. Untracked tutorial outputs are not copied.
`nwgmm` has the same character count as `dgemm`, so fixed-form columns do not move.
It is a separate symbol and is not rewritten by NWChem's BLAS conversion tools.
The modified link line places the shared backend after module archives and adds
its runtime library search path. Normal external BLAS/LAPACK dependencies remain.

**ABI configuration matters:** `NWGPU_FORTRAN_INTEGER_BYTES` must match the
NWChem compiler's default `INTEGER` and `LOGICAL` width (normally 8), independently
of `BLAS_SIZE`. The backend's CBLAS provider must use **32-bit BLAS integers**.
With an external LP64 BLAS, configure the normal NWChem `BLAS_SIZE=4` /
`USE_64TO32` path as appropriate. The bridge is tested with GNU Fortran's trailing
underscore naming and 4- and 8-byte defaults; other compiler conventions require
validation. Individual GEMM dimensions and leading dimensions are limited to
`INT_MAX`; oversized requests fail instead of narrowing silently.

## Run and measure

The full NWChem input format and executable come from `upstream/`. `nwgpu_bench`
is a kernel benchmark, not a molecular simulation executable.

| Variable | Meaning | Default |
| --- | --- | --- |
| `NWGPU_BACKEND` | `auto`, `cpu`, or `cuda` for Fortran calls | `auto` |
| `NWGPU_DEVICE` | CUDA device ordinal within `CUDA_VISIBLE_DEVICES`; set before first use | `0` |
| `NWGPU_MIN_GEMM_FLOPS` | Automatic GPU threshold, `2*m*n*k` | `100000000` |
| `NWGPU_MIN_GRID_POINTS` | Automatic exchange GPU threshold | `32768` |
| `NWGPU_REPORT_STATS` | Print per-process dispatch/fallback counts to stderr at normal exit when `1` | unset |

Automatic mode uses CPU for small operations or unavailable devices; a device
failure falls back to CPU before host outputs are changed. Explicit `cuda` fails
on an unavailable backend or device error. Trivial/empty GEMMs still use quick
paths. The public API returns a status; unrecoverable Fortran bridge errors print
a diagnostic and abort, rather than allowing invalid scientific results to continue.

For multiple GPUs, run the wrapper **inside** the MPI launcher so each local rank
selects a device. For example, with two GPUs and two local ranks:

```sh
mpirun -np 2 python3 tools/run_nwchem.py --gpus-per-node 2 --backend auto \
  -- /absolute/path/to/ported/nwchem input.nw
```

If the scheduler exposes only one GPU per rank through `CUDA_VISIBLE_DEVICES`,
use `--gpus-per-node 1`. Buffers and cuBLAS handles are cached per host thread;
multiple MPI ranks sharing a device consume separate device memory. No distributed
GPU matrix algorithm or automatic rank affinity inference is implemented.

Benchmark output reports **actual CPU/GPU dispatch counts**, elapsed time including
transfers, and error against CPU output. Compare matched CPU/GPU runs with the
same method, basis, grid, convergence settings, geometry and process count. First
compare energies, gradients, optimization convergence and frequencies against
original NWChem QA. Then measure complete wall time, peak host memory and device
memory. Kernel timings alone do not establish an application speedup. Neither
translation nor GPU use guarantees lower runtime or lower resource consumption.

## Remaining work for every module

See [MODULE_STATUS.md](MODULE_STATUS.md). Full migration requires validated C++
ports of each numerical subsystem and its supporting runtime, followed by the
appropriate module QA and cross-module workflow tests. This is a continuing
scientific software project; the current deliverable establishes an integrated
starting point while keeping all original modules available in source.

Ported NWChem routines retain the original ECL-2.0 licensing; see [LICENSE.TXT](LICENSE.TXT).
