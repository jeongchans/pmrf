#ifndef _NMF_H_
#define _NMF_H_

#include "common/numeric.h"

struct NMFMatrix {
    Float2dArray weight;
    Float2dArray basis;
};

class NMF {
  public:
    NMF(const double& tol=1e-10, const double& timelimit=0, const size_t& maxiter=8000, const int& verbose=0)
    : tol(tol), 
      timelimit(timelimit), 
      maxiter(maxiter), 
      verbose(verbose) {};

    void factorize(const Float2dArray& V, const int& r, Float2dArray& W, Float2dArray& H) const;
    Float1dArray transform(const Float1dArray& x, const Float2dArray& W) const;
    Float2dArray transform(const Float2dArray& V, const Float2dArray& W) const;

  private:
    double tol;
    double timelimit;
    size_t maxiter;
    int verbose;

    void nlssubprob(const Float2dArray V, const Float2dArray W,
                    const Float2dArray Hinit,
                    const double tol, const size_t maxiter,
                    Float2dArray& H, Float2dArray& grad, size_t& it) const;
    void nmf(const Float2dArray& V, const int& r,
             const Float2dArray& Winit, const Float2dArray& Hinit,
             Float2dArray& W, Float2dArray& H) const;
    Float2dArray nmf_feature(const Float2dArray& V, const Float2dArray& W) const;
};

///*
// * Projected gradient discriminant NMF
// */
//class PGDNMF {
//  public:
//    PGDNMF() : tol(1e-10), timelimit(120), maxiter(8000), verbose(0), gamma(1.0), delta(1.0) {};
//
//    NMFMatrix factorize(const Float2dArray& X, const int& M, const int& K, const Int1dArray& cls, const Float1dArray& wt);
//    Float1dArray transform(const Float1dArray& x, const Float2dArray& W) const;
//    Float2dArray transform(const Float2dArray& V, const Float2dArray& W) const;
//
//  private:
//    double tol;
//    double timelimit;
//    size_t maxiter;
//    int verbose;
//
//    double gamma;
//    double delta;
//
//    struct GradScatterMatrix {
//        Float2dArray grad_sw;
//        Float2dArray gradd_sw;
//        Float2dArray grad_sb;
//        Float2dArray gradd_sb;
//    };
//
//    struct SubSolution {
//        Float2dArray matrix;
//        Float2dArray grad;
//        size_t iter;
//    };
//
//    NMFMatrix optimize_projgrad(const Float2dArray& X, const Float2dArray& Zinit, const Float2dArray& Hinit, const int& K, const Int1dArray& cls);
//    GradScatterMatrix gradient_scatter(const Float2dArray& X, const Float2dArray& Z, const int& K, const Int1dArray& cls);
//    SubSolution solve_Z(const Float2dArray& X, const Float2dArray& Z, const Float2dArray& H, const Int1dArray& cls, const double& tol, const size_t& maxiter);
//    SubSolution solve_H(const Float2dArray& X, const Float2dArray& Z, const Float2dArray& H, const double& tol, const size_t& maxiter);
//};

#endif
