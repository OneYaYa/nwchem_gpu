# Migration status for every original source directory

This inventory is generated from `migration_manifest.json`. All source directories are retained in `upstream/src`; `(build)` groups files directly under `src/`.
Only DGEMM and xc_dirac/d2/d3 have C++/CUDA implementations. All other Fortran remains.
Counts are static call sites, including optional modules and generated sources; they do not establish runtime coverage.

| Directory | Tracked files | Fortran/include files | Routed GEMM sites | Current migration |
| --- | ---: | ---: | ---: | --- |
| `(build)` | 8 | 2 | 0 | Original implementation retained |
| `64to32blas` | 80 | 79 | 1 | GEMM routed; remaining Fortran retained |
| `NWints` | 448 | 412 | 32 | GEMM routed; remaining Fortran retained |
| `analyz` | 44 | 42 | 0 | Original implementation retained |
| `argos` | 759 | 758 | 0 | Original implementation retained |
| `atomscf` | 33 | 30 | 0 | Original implementation retained |
| `basis` | 1998 | 27 | 0 | Original implementation retained |
| `blas` | 113 | 109 | 0 | Original implementation retained |
| `bq` | 10 | 9 | 0 | Original implementation retained |
| `bsemol` | 27 | 26 | 14 | GEMM routed; remaining Fortran retained |
| `cafe` | 29 | 28 | 0 | Original implementation retained |
| `cca` | 43 | 1 | 0 | Original implementation retained |
| `ccca` | 3 | 2 | 0 | Original implementation retained |
| `cckohn` | 20 | 18 | 0 | Original implementation retained |
| `ccsd` | 70 | 63 | 139 | GEMM routed; remaining Fortran retained |
| `chelp` | 8 | 5 | 0 | Original implementation retained |
| `config` | 65 | 1 | 0 | Original implementation retained |
| `cons` | 9 | 8 | 0 | Original implementation retained |
| `cphf` | 10 | 8 | 0 | Original implementation retained |
| `dangchang` | 5 | 4 | 0 | Original implementation retained |
| `data` | 453 | 0 | 0 | Original implementation retained |
| `ddscf` | 160 | 154 | 0 | Original implementation retained |
| `develop` | 30 | 28 | 18 | GEMM routed; remaining Fortran retained |
| `diana` | 53 | 52 | 0 | Original implementation retained |
| `dimqm` | 54 | 53 | 1 | GEMM routed; remaining Fortran retained |
| `dntmc` | 5 | 4 | 0 | Original implementation retained |
| `dplot` | 15 | 14 | 0 | Original implementation retained |
| `drdy` | 61 | 59 | 0 | Original implementation retained |
| `driver` | 14 | 11 | 8 | GEMM routed; remaining Fortran retained |
| `esp` | 25 | 23 | 0 | Original implementation retained |
| `etrans` | 12 | 11 | 0 | Original implementation retained |
| `fcidump` | 4 | 3 | 0 | Original implementation retained |
| `fft` | 115 | 100 | 0 | Original implementation retained |
| `geninterface` | 6 | 3 | 0 | Original implementation retained |
| `geom` | 23 | 19 | 12 | GEMM routed; remaining Fortran retained |
| `gradients` | 11 | 10 | 0 | Original implementation retained |
| `gwmol` | 38 | 36 | 15 | GEMM routed; remaining Fortran retained |
| `hessian` | 42 | 39 | 7 | GEMM routed; remaining Fortran retained |
| `include` | 2 | 0 | 0 | Original implementation retained |
| `inp` | 13 | 8 | 0 | Original implementation retained |
| `input` | 12 | 9 | 0 | Original implementation retained |
| `java` | 58 | 0 | 0 | Original implementation retained |
| `lapack` | 464 | 461 | 118 | GEMM routed; remaining Fortran retained |
| `leps` | 23 | 19 | 0 | Original implementation retained |
| `libext` | 40 | 1 | 0 | Original implementation retained |
| `lucia` | 166 | 81 | 47 | GEMM routed; remaining Fortran retained |
| `mcscf` | 51 | 47 | 27 | GEMM routed; remaining Fortran retained |
| `mm` | 29 | 23 | 0 | Original implementation retained |
| `moints` | 26 | 22 | 25 | GEMM routed; remaining Fortran retained |
| `mp2_grad` | 24 | 22 | 43 | GEMM routed; remaining Fortran retained |
| `mrpt` | 76 | 58 | 3 | GEMM routed; remaining Fortran retained |
| `nbo` | 1 | 0 | 0 | Original implementation retained |
| `ncc` | 8 | 7 | 0 | Original implementation retained |
| `nwc_columbus` | 36 | 30 | 0 | Original implementation retained |
| `nwdft` | 621 | 589 | 35 | Dirac d1/d2/d3 ported; GEMM routed; remaining Fortran retained |
| `nwmd` | 11 | 10 | 0 | Original implementation retained |
| `nwpw` | 1687 | 957 | 76 | GEMM routed; remaining Fortran retained |
| `nwxc` | 511 | 299 | 0 | Original implementation retained |
| `optim` | 25 | 19 | 12 | GEMM routed; remaining Fortran retained |
| `peigs` | 212 | 36 | 0 | Original implementation retained |
| `peigs_comm` | 12 | 6 | 0 | Original implementation retained |
| `perfm` | 10 | 7 | 0 | Original implementation retained |
| `prepar` | 112 | 110 | 0 | Original implementation retained |
| `property` | 80 | 79 | 0 | Original implementation retained |
| `pstat` | 13 | 12 | 0 | Original implementation retained |
| `python` | 14 | 1 | 0 | Original implementation retained |
| `qhop` | 13 | 11 | 0 | Original implementation retained |
| `qmd` | 44 | 41 | 0 | Original implementation retained |
| `qmmm` | 34 | 33 | 0 | Original implementation retained |
| `rdmft` | 293 | 270 | 199 | GEMM routed; remaining Fortran retained |
| `rimp2` | 38 | 34 | 2 | GEMM routed; remaining Fortran retained |
| `rimp2_grad` | 55 | 37 | 2 | GEMM routed; remaining Fortran retained |
| `rism` | 74 | 11 | 1 | GEMM routed; remaining Fortran retained |
| `rtdb` | 31 | 6 | 0 | Original implementation retained |
| `selci` | 101 | 96 | 2 | GEMM routed; remaining Fortran retained |
| `smd` | 591 | 456 | 0 | Original implementation retained |
| `solvation` | 22 | 21 | 0 | Original implementation retained |
| `space` | 4 | 3 | 0 | Original implementation retained |
| `stepper` | 53 | 47 | 1 | GEMM routed; remaining Fortran retained |
| `symmetry` | 50 | 47 | 5 | GEMM routed; remaining Fortran retained |
| `task` | 27 | 25 | 1 | GEMM routed; remaining Fortran retained |
| `tce` | 11941 | 11648 | 12236 | GEMM routed; remaining Fortran retained |
| `tools` | 17 | 0 | 0 | Original implementation retained |
| `uccsdt` | 26 | 25 | 39 | GEMM routed; remaining Fortran retained |
| `util` | 254 | 172 | 4 | GEMM routed; remaining Fortran retained |
| `vib` | 33 | 32 | 15 | GEMM routed; remaining Fortran retained |
| `vscf` | 3 | 2 | 0 | Original implementation retained |
| `xtb` | 11 | 10 | 0 | Original implementation retained |

## Work needed to reach a complete C++/CUDA application

1. Establish reproducible full CPU builds and module QA for the original and modified trees. Capture method, basis, task, process count, energies, derivatives, convergence, time and memory.
2. Profile representative workloads for every enabled module. Use these measurements to choose subsequent ports and GPU batching boundaries; a static GEMM count is not a profile.
3. Port and validate Gaussian basis evaluation, one-/two-electron integrals, screening, Coulomb and exact-exchange construction. Preserve symmetry and integral conventions.
4. Port the remaining exchange-correlation functionals and their response/gradient kernels, then density/grid construction. Keep intermediate grids and matrices resident on device across calls.
5. Port tensor contraction/permutation operations and distributed data movement for correlated methods. Account for the existing TCE CUDA implementation rather than assuming all TCE work is CPU-only.
6. Port FFT/plane-wave operations and validate pseudopotential, PAW, periodic, spin, and complex arithmetic paths. Current DGEMM routing does not accelerate complex GEMM or FFTs.
7. Port molecular dynamics forces, gradients, Hessians, response, optimization and vibrations, validating complete workflows as well as isolated kernels.
8. Replace the remaining Fortran runtime, parsers, scientific I/O and supporting modules with C++ as independently reviewed migrations. Preserve input compatibility and scientific data formats.
9. Validate MPI/Global Arrays behavior, rank/device mapping, memory use, failure behavior and multi-GPU scaling. Complete original module QA and document unsupported configurations before claiming feature parity.

No stage beyond the current kernel port/source routing is recorded as complete. No GPU execution or full-application validation has been performed on this Mac.
