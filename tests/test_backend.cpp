#include "nwchem_gpu.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void expect(bool ok, const char* message) { if (!ok) throw std::runtime_error(message); }
void close(double actual, double expected, double tolerance = 2e-12) {
    if (!(std::abs(actual-expected) <= tolerance * std::max(1., std::abs(expected))))
        throw std::runtime_error("numerical mismatch: " + std::to_string(actual) + " vs " + std::to_string(expected));
}
void success(int result) { if (result) throw std::runtime_error(nwgpu_last_error()); }
bool trans(char op) { return op != 'N'; }
void matrix_tests(int backend) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> random(-1., 1.);
    for (char ta : {'N', 'T', 'C'}) for (char tb : {'N', 'T', 'C'})
    for (double alpha : {1.3, -0.75}) for (double beta : {0., -.25, 1.}) {
        constexpr int m=7, n=5, k=9, ldc=m+3;
        const int ar=trans(ta)?k:m, ac=trans(ta)?m:k;
        const int br=trans(tb)?n:k, bc=trans(tb)?k:n;
        const int lda=ar+2, ldb=br+4;
        std::vector<double> a(lda*ac), b(ldb*bc), c(ldc*n, 12345.), expected=c;
        for (double& v:a) v=random(rng);
        for (double& v:b) v=random(rng);
        for (int j=0;j<n;++j) for (int i=0;i<m;++i) {
            const int idx=j*ldc+i;
            const double initial=random(rng);
            // NaN detects an erroneous read of C with beta=0.
            c[idx]=beta==0 ? std::numeric_limits<double>::quiet_NaN() : initial;
            double sum=0;
            for (int p=0;p<k;++p)
                sum += a[trans(ta)?i*lda+p:p*lda+i] * b[trans(tb)?p*ldb+j:j*ldb+p];
            expected[idx]=alpha*sum+(beta==0?0:beta*initial);
        }
        success(nwgpu_dgemm(backend,ta,tb,m,n,k,alpha,a.data(),lda,b.data(),ldb,beta,c.data(),ldc));
        for (size_t i=0;i<c.size();++i) close(c[i],expected[i]);
    }
    double c=std::numeric_limits<double>::quiet_NaN();
    success(nwgpu_dgemm(backend,'N','N',1,1,3,0,nullptr,1,nullptr,3,0,&c,1)); close(c,0);
    c=3;
    success(nwgpu_dgemm(backend,'N','N',1,1,0,1,nullptr,1,nullptr,1,2,&c,1)); close(c,6);
    success(nwgpu_dgemm(backend,'N','N',0,4,3,1,nullptr,1,nullptr,3,0,nullptr,1));
    expect(nwgpu_dgemm(backend,'X','N',1,1,1,1,&c,1,&c,1,0,&c,1)!=0,"transpose validation");
    expect(nwgpu_dgemm(backend,'N','N',2,1,1,1,&c,1,&c,1,0,&c,2)!=0,"leading dimension validation");
    expect(nwgpu_dgemm(backend,'N','N',INT64_MAX,1,1,1,&c,1,&c,1,0,&c,1)!=0,"integer narrowing validation");
    expect(nwgpu_dgemm(backend,'N','N',INT32_MAX,INT32_MAX,1,1,&c,INT32_MAX,&c,1,0,&c,INT32_MAX)!=0,"matrix extent overflow validation");
    expect(nwgpu_dgemm(99,'N','N',1,1,1,1,&c,1,&c,1,0,&c,1)!=0,"backend validation");
}
void dirac_tests(int backend) {
    // Analytic r=1 and r=8, strict cutoff, negative density, seeded outputs.
    const double tol=1e-12, fac=.7, p1=-.9847450218426959, p2=-.7385587663820219;
    constexpr int n=5;
    double rho[n]={1,8,0,-1,tol}, weights[n]={.1,.2,.3,.4,.5};
    std::vector<double> a(n,.2), a2(3*n,.3), a3(4*n,.4), f(n,.5);
    double ex=2.;
    success(nwgpu_dirac(backend,3,n,1,tol,fac,rho,weights,a.data(),a2.data(),a3.data(),&ex,1,f.data()));
    close(ex,2.+p2*fac*(.1+16*.2));
    for (int i=0;i<n;++i) {
        const double root=i==0?1:i==1?2:0;
        close(a[i],(.2+root)*p1*fac);
        close(a2[i],(.3+(i==0?1:i==1?.25:0))*(2./3.)*p1*fac);
        close(a3[i],(.4+(i==0?1:i==1?1./32.:0))*(-8./9.)*p1*fac);
        close(f[i],.5+(i<2?rho[i]*root*p2*fac:0));
    }
    for (int i=n;i<3*n;++i) close(a2[i],.3);
    for (int i=n;i<4*n;++i) close(a3[i],.4);
    // Unrestricted convention: total density then alpha and beta, one empty spin.
    const double srho[6]={9,1,1,1,8,0}, w[2]={.25,.75};
    double sa[4]={0}, s2[6]={0}, s3[8]={0}, sf[2]={0};
    ex=0;
    success(nwgpu_dirac(backend,3,2,2,tol,1,srho,w,sa,s2,s3,&ex,1,sf));
    const double p3=-1.240700981798799, p4=-.9305257363490993;
    close(ex,p4*(17*.25+1*.75));
    close(sa[0],p3); close(sa[1],p3); close(sa[2],2*p3); close(sa[3],0);
    close(s2[0],p3/3); close(s2[4],p3/12); close(s2[5],0);
    close(s3[0],-2*p3/9); close(s3[6],-p3/144); close(s3[7],0);
    // Potential is the derivative of the energy; kernel is derivative of potential.
    for (int spin=0;spin<2;++spin) {
        const double r=.8+spin, h=1e-5;
        double energy[2], potential[2];
        for (int sign=0;sign<2;++sign) {
            const double v=r+(sign?1:-1)*h;
            const double rr[3]={v+.6,spin==0?v:.6,spin==1?v:.6}, ww=1;
            double aa[2]={0}, ee=0;
            success(nwgpu_dirac(backend,1,1,2,tol,1,rr,&ww,aa,nullptr,nullptr,&ee,0,nullptr));
            energy[sign]=ee; potential[sign]=aa[spin];
        }
        close((energy[1]-energy[0])/(2*h),p3*std::cbrt(r),1e-9);
        close((potential[1]-potential[0])/(2*h),p3/(3*std::pow(r,2./3.)),1e-9);
    }
    // A block boundary and a large grid exercise GPU reduction/tail coverage.
    for (int q : {1,255,256,257,33001}) {
        std::vector<double> r(q,1.), ww(q,1./q), aa(q,0.);
        double ee=0;
        success(nwgpu_dirac(backend,1,q,1,tol,1,r.data(),ww.data(),aa.data(),nullptr,nullptr,&ee,0,nullptr));
        close(ee,p2); for (double v:aa) close(v,p1);
    }
    expect(nwgpu_dirac(backend,2,1,1,tol,1,rho,weights,a.data(),nullptr,nullptr,&ex,0,nullptr)!=0,"derivative storage validation");
}
}
#if NWGPU_FORTRAN_INTEGER_BYTES == 8
using FI=int64_t;
#else
using FI=int32_t;
#endif
extern "C" void nwgmm_(const char*,const char*,const FI*,const FI*,const FI*,const double*,
                        const double*,const FI*,const double*,const FI*,const double*,double*,const FI*);
int main(int argc,char** argv) {
    try {
        const int backend=argc>1 && std::string(argv[1])=="cuda" ? NWGPU_CUDA : NWGPU_CPU;
        matrix_tests(backend); dirac_tests(backend);
        const FI one=1; const char nn='N'; const double a=2,b=3,alpha=1,beta=0; double c=0;
        nwgmm_(&nn,&nn,&one,&one,&one,&alpha,&a,&one,&b,&one,&beta,&c,&one); close(c,6);
        if (backend==NWGPU_CUDA) {
            const auto stats=nwgpu_get_stats();
            expect(stats.cuda_gemm>=54 && stats.cuda_dirac>=7,"CUDA tests did not run on device");
        } else if (!nwgpu_cuda_available()) {
            expect(nwgpu_dgemm(NWGPU_CUDA,'N','N',1,1,1,1,&a,1,&b,1,0,&c,1)!=0,"strict CUDA must fail without device");
        }
        std::cout<<"All GEMM, Dirac, derivative, ABI and error contracts passed\n";
    } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
