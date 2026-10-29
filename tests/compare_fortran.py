#!/usr/bin/env python3
"""Compile the original NWChem routines, then compare through their real ABI.

No NumPy needed. Requires gfortran; creates all oracle artifacts under --work.
"""
import argparse
import ctypes as ct
import json
import math
import os
from pathlib import Path
import random
import subprocess


def main():
    project = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", required=True)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--source", type=Path, default=project.parent)
    parser.add_argument("--work", type=Path, default=project / "build-oracle")
    parser.add_argument("--integer-bytes", type=int, choices=(4, 8), default=8)
    parser.add_argument("--backend", choices=("cpu", "cuda"), default="cpu")
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    (work / "errquit.fh").write_text("! No error calls in the standalone Dirac oracle.\n")
    (work / "stdio.fh").write_text("! No I/O in the standalone Dirac oracle.\n")
    (work / "support.f90").write_text("""subroutine dscal(n,alpha,x,incx)
implicit none
integer n,incx,i
double precision alpha,x(*)
do i=1,n
  x(1+(i-1)*incx)=alpha*x(1+(i-1)*incx)
enddo
end subroutine
logical function lsame(a,b)
implicit none
character a,b
integer ia,ib
ia=iachar(a)
ib=iachar(b)
if(ia>=97.and.ia<=122) ia=ia-32
if(ib>=97.and.ib<=122) ib=ib-32
lsame=ia==ib
end function
subroutine xerbla(name,info)
implicit none
character*(*) name
integer info
print *, 'Oracle invalid call', name, info
stop 1
end subroutine
""")
    suffix = ".dylib" if __import__("sys").platform == "darwin" else ".so"
    oracle_path = work / ("liboriginal" + suffix)
    command = [args.fc, "-shared", "-fPIC", "-cpp", "-O2", "-I" + str(work),
               "-I" + str(args.source.resolve() / "src/nwdft/include")]
    if args.integer_bytes == 8:
        command += ["-fdefault-integer-8"]
    command += [str(args.source.resolve() / "src/nwdft/xc/xc_dirac.F"),
                str(args.source.resolve() / "src/blas/double/dgemm.f"),
                str(work / "support.f90"), "-o", str(oracle_path)]
    subprocess.run(command, check=True)
    original = ct.CDLL(str(oracle_path))
    port = ct.CDLL(str(args.library.resolve()))
    fi = ct.c_int64 if args.integer_bytes == 8 else ct.c_int32
    ptr = ct.POINTER(ct.c_double)
    port.nwgpu_dirac.argtypes = [ct.c_int, ct.c_int, ct.c_int64, ct.c_int, ct.c_double,
                                ct.c_double, ptr, ptr, ptr, ptr, ptr, ptr, ct.c_int, ptr]
    port.nwgpu_last_error.restype = ct.c_char_p
    backend = 1 if args.backend == "cpu" else 2
    os.environ["NWGPU_BACKEND"] = args.backend
    rng = random.Random(20261029)
    worst = 0.
    cases = 0

    def array(values):
        return (ct.c_double * len(values))(*values)

    def compare(actual, expected):
        nonlocal worst
        for a, b in zip(actual, expected):
            error = abs(a-b) / max(1., abs(b))
            worst = max(worst, error)
            if not math.isfinite(error) or error > 5e-12:
                raise AssertionError((a, b, error))

    for ipol in (1, 2):
        for order in (1, 2, 3):
            for nq in (1, 17, 257):
                for ldew in (0, 1):
                    for fac in (0., .73, -1.):
                        tol = 1e-12
                        densities = [10**rng.uniform(-8, 2) for _ in range(nq*(1 if ipol == 1 else 3))]
                        for i in range(min(nq, 5)):
                            densities[i] = [0., -1., tol, tol*.5, tol*2][i]
                        if ipol == 2:
                            densities[2*nq] = 0. # fully polarized point
                        rho = array(densities)
                        weights = array([rng.uniform(.001, 1) for _ in range(nq)])
                        initial = [[rng.uniform(-.2, .2) for _ in range(nq*c)] for c in (ipol, 3, 4, 1)]
                        output = [array(x) for x in initial]
                        expected = [array(x) for x in initial]
                        ex, ex_ref = ct.c_double(.42), ct.c_double(.42)
                        f_tol, f_fac, f_q, f_pol, logical = ct.c_double(tol), ct.c_double(fac), fi(nq), fi(ipol), fi(ldew)
                        ignored = fi(0)
                        routine = getattr(original, "xc_dirac" + ("" if order == 1 else "_d" + str(order)) + "_")
                        call = [ct.byref(f_tol), ct.byref(f_fac), ct.byref(ignored), ct.byref(ignored), rho, expected[0]]
                        if order >= 2: call.append(expected[1])
                        if order == 3: call.append(expected[2])
                        call += [ct.byref(f_q), ct.byref(f_pol), ct.byref(ex_ref), weights, ct.byref(logical), expected[3]]
                        routine(*call)
                        result = port.nwgpu_dirac(backend, order, nq, ipol, tol, fac, rho, weights,
                                                   output[0], output[1], output[2], ct.byref(ex), ldew, output[3])
                        if result: raise RuntimeError(port.nwgpu_last_error().decode())
                        for actual, reference in zip(output, expected): compare(actual, reference)
                        compare([ex.value], [ex_ref.value])
                        # Exercise the Fortran-facing C++ entry point with identical ABI.
                        bridge = getattr(port, "xc_dirac" + ("" if order == 1 else "_d" + str(order)) + "_")
                        bridge_output = [array(x) for x in initial]
                        bridge_ex = ct.c_double(.42)
                        bridge_call = call[:5] + [bridge_output[0]]
                        if order >= 2: bridge_call.append(bridge_output[1])
                        if order == 3: bridge_call.append(bridge_output[2])
                        bridge_call += [ct.byref(f_q), ct.byref(f_pol), ct.byref(bridge_ex), weights, ct.byref(logical), bridge_output[3]]
                        bridge(*bridge_call)
                        for actual, reference in zip(bridge_output, expected): compare(actual, reference)
                        compare([bridge_ex.value], [ex_ref.value])
                        cases += 1
    # Compare the routed Fortran GEMM entry point to the original blocked DGEMM.
    for ta in (b"N", b"T", b"C"):
        for tb in (b"N", b"T", b"C"):
            for beta in (0., -.4, 1.):
                m, n, k = 7, 5, 9
                lda, ldb, ldc = (m if ta == b"N" else k)+2, (k if tb == b"N" else n)+3, m+4
                a = array([rng.random() for _ in range(lda*(k if ta == b"N" else m))])
                b = array([rng.random() for _ in range(ldb*(n if tb == b"N" else k))])
                c = array([rng.random() for _ in range(ldc*n)])
                reference = array(list(c))
                ints = [fi(x) for x in (m,n,k,lda,ldb,ldc)]
                alpha, fbeta = ct.c_double(.73), ct.c_double(beta)
                call = [ct.c_char_p(ta),ct.c_char_p(tb),*[ct.byref(x) for x in ints[:3]],
                        ct.byref(alpha),a,ct.byref(ints[3]),b,ct.byref(ints[4]),ct.byref(fbeta)]
                # gfortran CHARACTER hidden lengths are trailing size_t arguments.
                original.dgemm_(*call,reference,ct.byref(ints[5]),ct.c_size_t(1),ct.c_size_t(1))
                port.nwgmm_(*call,c,ct.byref(ints[5]),ct.c_size_t(1),ct.c_size_t(1))
                compare(c,reference)
                cases += 1
    result = {"backend": args.backend, "integer_bytes": args.integer_bytes,
              "original_fortran_comparison_cases": cases, "max_scaled_error": worst,
              "scope": "Dirac d1/d2/d3 and DGEMM, including C++ Fortran ABI"}
    (work / "result.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
