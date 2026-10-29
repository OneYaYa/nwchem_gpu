#pragma once
#include "internal.hpp"
#ifdef __CUDACC__
#define NWGPU_HD __host__ __device__
#else
#define NWGPU_HD
#endif
namespace nwgpu {
// Translated from src/nwdft/xc/xc_dirac.F. Original ECL-2.0 notice in LICENSE.TXT.
// Each point owns its output locations. Return the unscaled energy contribution.
NWGPU_HD inline double dirac_point(const Dirac& d, int i) {
    constexpr double p1 = -0.9847450218426959;
    constexpr double p2 = -0.7385587663820219;
    constexpr double p3 = -1.240700981798799;
    constexpr double p4 = -0.9305257363490993;
    const size_t q = static_cast<size_t>(d.nq);
    double term = 0.;
    if (d.ipol == 1) {
        const double r = d.rho[i];
        if (r > d.tol) {
            const double root = pow(r, 1. / 3.);
            d.a[i] += root;
            if (d.order >= 2) d.a2[i] += 1. / pow(r, 2. / 3.);
            if (d.order == 3) d.a3[i] += 1. / pow(r, 5. / 3.);
            term = r * root * d.w[i];
            if (d.ldew) d.func[i] += r * root * d.fac * p2;
        }
        // The Fortran routine scales even pre-existing/screened Amat entries.
        d.a[i] *= p1 * d.fac;
        if (d.order >= 2) d.a2[i] *= (2. / 3.) * p1 * d.fac;
        if (d.order == 3) d.a3[i] *= (-8. / 9.) * p1 * d.fac;
    } else if (d.rho[i] > d.tol) {
        for (size_t spin = 0; spin < 2; ++spin) {
            const double r = d.rho[(spin + 1) * q + i];
            const double root = r > d.tol ? pow(r, 1. / 3.) : 0.;
            d.a[spin * q + i] += p3 * root * d.fac;
            if (r > d.tol) {
                if (d.order >= 2)
                    d.a2[2 * spin * q + i] += (1. / 3.) * p3 * d.fac / pow(r, 2. / 3.);
                if (d.order == 3)
                    d.a3[3 * spin * q + i] -= (2. / 9.) * p3 * root / (r * r) * d.fac;
            }
            term += root * r;
        }
        if (d.ldew) d.func[i] += term * p4 * d.fac;
        term *= d.w[i];
    }
    return term;
}
inline double dirac_energy_factor(const Dirac& d) {
    return (d.ipol == 1 ? -0.7385587663820219 : -0.9305257363490993) * d.fac;
}
}
#undef NWGPU_HD
