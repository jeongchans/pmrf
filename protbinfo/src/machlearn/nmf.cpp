/* NMF by alternative non-negative least squares using projected gradients
 * Matlab codes implemented by Chih-Jen Lin, National Taiwan Univ.
 * Ported to C++ by Chan-Seok Jeong
 */

#include "nmf.h"

#include <algorithm>
#include <iostream>
#include <ctime>

#include <gsl/gsl_linalg.h>

using blitz::pow2;
using blitz::sum;
using blitz::where;

inline Float2dArray _omp_dot(const Float2dArray& x, const Float2dArray& y) {
    int n = x.rows();
    int l = x.cols();
    int m = y.cols();
    Float2dArray r = zeros(n, m);
    int i, j, k;
#pragma omp parallel for private(j, k) shared(r, x, y)
    for (i = 0; i < n; ++i) {
        for (j = 0; j < m; ++j) {
            for (k = 0; k < l; ++k) r(i, j) += x(i, k) * y(k, j);
        }
    }
    return r;
}

inline Float2dArray _dot(const Float2dArray& x, const Float2dArray& y) {
    return _omp_dot(x, y);
}

Float2dArray _pinv_matrix(const Float2dArray& A) {
    Float2dArray A_dup(A.shape());
    A_dup = A;
    int rows = A.rows();
    int cols = A.cols();
    gsl_matrix_view m = gsl_matrix_view_array(A_dup.data(), rows, cols);
    gsl_matrix *v = gsl_matrix_alloc(cols, cols);
    gsl_vector *s = gsl_vector_alloc(cols);
    gsl_linalg_SV_decomp_jacobi(&m.matrix, v, s);
    Float2dArray U(rows, cols);
    Float2dArray S(cols, cols);
    Float2dArray V(cols, cols);
    Float2dArray St(cols, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j) U(i, j) = gsl_matrix_get(&m.matrix, i, j);
    for (int i = 0; i < cols; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (i == j) S(i, j) = gsl_vector_get(s, j);
            else S(i, j) = 0;
            if (S(i, j) != 0) St(j, i) = 1. / S(i, j);
            else St(i, j) = 0;
        }
    }
    for (int i = 0; i < cols; ++i)
        for (int j = 0; j < cols; ++j) V(i, j) = gsl_matrix_get(v, i, j);
    Float2dArray r(cols, rows);
    r = dot(dot(V, St), transpose(U));
    gsl_matrix_free(v);
    gsl_vector_free(s);
    return r;
}

void NMF::factorize(const Float2dArray& V, const int& r, Float2dArray& W, Float2dArray& H) const {
    int m = V.rows();
    int n = V.cols();
    Float2dArray Winit(n, r);
    Winit = blitz::abs(randn(n, r));
    Float2dArray Hinit(r, m);
    Hinit = blitz::abs(randn(r, m));
    nmf(transpose(V), r, Winit, Hinit, W, H);
    H.transposeSelf(blitz::secondDim, blitz::firstDim);
}

Float1dArray NMF::transform(const Float1dArray& x, const Float2dArray& W) const {
    Float2dArray V(1, x.size());
    V(0, ALL) = x.copy();
    return transform(V, W)(0, ALL);
}

Float2dArray NMF::transform(const Float2dArray& V, const Float2dArray& W) const {
    Float2dArray H = nmf_feature(transpose(V), W);
    H.transposeSelf(blitz::secondDim, blitz::firstDim);
    return H;
}

/* NMF by alternative non-negative least squares using projected gradients
 *
 * Args:    V: target matrix
 *          r: number of basis
 *          Winit: initial solution
 *          Hinit: initial solution
 *
 * Returns: W: solution
 *          H: solution
 */

void NMF::nmf(const Float2dArray& V, const int& r,
              const Float2dArray& Winit, const Float2dArray& Hinit,
              Float2dArray& W, Float2dArray& H) const {
    double projnorm;
    size_t it, itW, itH;
    time_t initt = time(NULL);
    int n = V.rows();
    int m = V.cols();
    W.resize(Winit.shape());
    W = Winit.copy();
    H.resize(Hinit.shape());
    H = Hinit.copy();
    Float2dArray gradW(n, r);
    gradW = _dot(W, _dot(H, transpose(H))) - _dot(V, transpose(H));
    Float2dArray gradH(r, m);
    gradH = _dot(_dot(transpose(W), W), H) - _dot(transpose(W), V);
    double initgrad = norm(concatenate(gradW, transpose(gradH)));
    if (verbose > 0) std::cout << "Init gradient norm " << initgrad << std::endl;
    double tolW = std::max(0.001, tol) * initgrad;
    double tolH = tolW;
    for (it = 0; it < maxiter; ++it) {
        projnorm = sqrt(sum(where(gradW < 0 || W > 0, pow2(gradW), 0)) + 
                        sum(where(gradH < 0 || H > 0, pow2(gradH), 0)));
        if (projnorm < tol * initgrad) break;
        if (timelimit > 0 && difftime(time(NULL), initt) > timelimit) break;
        if (verbose > 0) {
            std::cout << "Iter = " << it << std::endl;
            std::cout << "Current proj-grad norm " << projnorm << std::endl;
        }
        nlssubprob(transpose(V), transpose(H), transpose(W), tolW, 1000, W, gradW, itW);
        W.transposeSelf(blitz::secondDim, blitz::firstDim);
        gradW.transposeSelf(blitz::secondDim, blitz::firstDim);
        if (itW == 0) tolW = 0.1 * tolW;
        nlssubprob(V, W, H, tolH, 1000, H, gradH, itH);
        if (itH == 0) tolH = 0.1 * tolH;
        if (verbose > 0 && it % 10 == 0) {
            std::cout << '.';
            std::cout.flush();
        }
    }
    if (verbose > 0) {
        std::cout << "Iter = " << it << std::endl;
        std::cout << "Final proj-grad norm " << projnorm << std::endl;
    }
}

Float2dArray NMF::nmf_feature(const Float2dArray& V, const Float2dArray& W) const {
    int r = W.cols();
    int m = V.cols();
    Float2dArray H(r, m);
    //H = blitz::abs(randn(r, m));
    H = dot(transpose(W), V);
    Float2dArray gradH(r, m);
    size_t itH;
    nlssubprob(V, W, H, tol, maxiter, H, gradH, itH);
    return H;
}

//Float2dArray NMF::nmf_feature(const Float2dArray& V, const Float2dArray& W) const {
//    double projnorm, prev_projnorm;
//    size_t it, itH;
//    int r = W.cols();
//    int n = V.rows();
//    int m = V.cols();
//    Float2dArray H(r, m);
//    //H = blitz::abs(randn(r, m));
//    H = dot(transpose(W), V);
//    //H = dot(_pinv_matrix(W), V);
//    //for (int i = 0; i < r; ++i) {
//    //    for (int j = 0; j < m; ++j) {
//    //        if (H(i, j) < 0) H(i, j) = 0;
//    //    }
//    //}
//    Float2dArray gradH(r, m);
//    gradH = dot(dot(transpose(W), W), H) - dot(transpose(W), V);
//    double initgrad = norm(transpose(gradH));
//    //double tolH = 0.001 * initgrad;
//    double tolH = tol * initgrad;
//    //nlssubprob(V, W, H, tolH, 1000, H, gradH, itH);
//    //nlssubprob(V, W, H, tolH, maxiter, H, gradH, itH);
//    nlssubprob(V, W, H, tol, maxiter, H, gradH, itH);
//    return H;
//
//    for (it = 0; it < maxiter; ++it) {
//        projnorm = sqrt(sum(where(gradH < 0 || H > 0, pow2(gradH), 0)));
//        if (projnorm < tol * initgrad) break;
//        nlssubprob(V, W, H, tolH, 1000, H, gradH, itH);
//        if (itH == 0) tolH = 0.1 * tolH;
//    }
//    return H;
//}

/* Solving the sub-problem by the projected gradient algorithm
 *
 * Args:    V, W: constant matrics
 *          Hinit: initial solution
 *          tol: stopping tolerance
 *          maxiter: limit of iterations
 *
 * Returns: H: resulted solution
 *          grad: gradient
 *          it: iteration
 */

void NMF::nlssubprob(const Float2dArray V, const Float2dArray W,
                     const Float2dArray Hinit,
                     const double tol, const size_t maxiter,
                     Float2dArray& H, Float2dArray& grad, size_t& it) const {
    Float2dArray Hn(Hinit.shape());
    Float2dArray Hp(Hinit.shape());
    Float2dArray d(Hinit.shape());
    double projgrad, gradd, dQd;
    bool suff_decr, decr_alpha;
    H.resize(Hinit.shape());
    H = Hinit.copy();
    grad.resize(Hinit.shape());
    Float2dArray WtV = _dot(transpose(W), V);
    Float2dArray WtW = _dot(transpose(W), W);
    double alpha = 1.0;
    double beta = 0.1;
    for (it = 0; it < maxiter; ++it) {
        grad = _dot(WtW, H) - WtV;
        //grad = dot(WtW, H) - WtV;
        projgrad = sqrt(sum(where(grad < 0 || H > 0, pow2(grad), 0)));
        if (projgrad < tol) break;
        for (int inner_it = 0; inner_it < 20; ++inner_it) {
            Hn = H - alpha * grad;
            Hn = where(Hn > 0, Hn.copy(), 0);
            d = Hn - H;
            gradd = sum(grad * d);
            dQd = sum(_dot(WtW, d) * d);
            //dQd = sum(dot(WtW, d) * d);
            suff_decr = (0.99 * gradd + 0.5 * dQd) < 0;
            if (inner_it == 0) {
                decr_alpha = !suff_decr;
                Hp = H.copy();
            }
            if (decr_alpha) {
                if (suff_decr) {
                    H = Hn.copy();
                    break;
                } else {
                    alpha = alpha * beta;
                }
            } else {
                if (!suff_decr || all(Hp == Hn)) {
                    H = Hp.copy();
                    break;
                } else {
                    alpha = alpha / beta;
                    Hp = Hn.copy();
                }
            }
        }
    }
    if (it >= maxiter) std::cerr << "Max iter in nlssubprob" << std::endl;
}

///*
// * Projected gradient discriminant NMF
// *
// * Refenrence:
// *   Lin. Projected gradient methods for nonnegative matrix factorization. 
// *   Neural Computation (2007) vol. 19 (10) pp. 2756-2779
// *
// * Args:    X: target matrix
// *          M: number of basis
// *          K: number of classes
// *          cls: classification[0, 1, 2, ...]
// *          wt: weight of class
// *
// * Returns: NMFMatrix
// */
//NMFMatrix PGDNMF::factorize(const Float2dArray& X, const int& M, const int& K, const Int1dArray& cls, const Float1dArray& wt) {
//    int F = X.rows();
//    int T = X.cols();
//    Float2dArray Zinit(F, M);  // initial weight matrix
//    Float2dArray Hinit(M, T);  // initial basis matrix
//    Zinit = blitz::abs(randn(F, M));
//    Hinit = blitz::abs(randn(M, T));
//    NMFMatrix m = optimize_projgrad(X, Zinit, Hinit, K, cls);
//    return m;
//}
//
//Float1dArray PGDNMF::transform(const Float1dArray& x, const Float2dArray& W) const {
//}
//
//Float2dArray PGDNMF::transform(const Float2dArray& V, const Float2dArray& W) const {
//}
//
//NMFMatrix PGDNMF::optimize_projgrad(const Float2dArray& X, const Float2dArray& Zinit, const Float2dArray& Hinit, const int& K, const Int1dArray& cls) {
//    Float2dArray H = Hinit.copy();
//    Float2dArray Z = Zinit.copy();
//    GradScatterMatrix gradS = gradient_scatter(X, Z, K, cls);
//    Float2dArray gradH(H.shape());
//    Float2dArray gradZ(Z.shape());
//    gradH = dot(transpose(Z), Float2dArray(dot(Z, H) - X));
//    gradZ = dot(Float2dArray(dot(Z, H) - X), transpose(H)) 
//            + gamma * gradS.grad_sw
//            - delta * gradS.grad_sb;
//    double initgrad = norm(gradH) + norm(gradZ);
//    if (verbose > 0) std::cout << "Init gradient norm " << initgrad << std::endl;
//    double tolZ = std::max(0.001, tol) * initgrad;
//    double tolH = tolZ;
//    int it, itH, itZ;
//    time_t initt = time(NULL);
//    double projnorm;
//    for (it = 0; it < maxiter; ++it) {
//        projnorm = sqrt(sum(where(gradZ < 0 || Z > 0, pow2(gradZ), 0)) + 
//                        sum(where(gradH < 0 || H > 0, pow2(gradH), 0)));
//        if (projnorm < tol * initgrad || difftime(time(NULL), initt) > timelimit)
//            break;
//        SubSolution sol1 = solve_Z(X, Z, H, cls, tolZ, 1000);
//        Z = sol1.matrix;
//        gradZ = sol1.grad;
//        itZ = sol1.iter;
//        if (itZ == 0) tolZ = 0.1 * tolZ;
//        SubSolution sol2 = solve_H(X, Z, H, tolH, 1000);
//        H = sol1.matrix;
//        gradH = sol1.grad;
//        itH = sol1.iter;
//        if (itH == 0) tolH = 0.1 * tolH;
//        if (verbose > 0 && it % 10 == 0) {
//            std::cout << '.';
//            std::cout.flush();
//        }
//    }
//    if (verbose > 0) {
//        std::cout << "Iter = " << it << std::endl;
//        std::cout << "Final proj-grad norm " << projnorm << std::endl;
//    }
//}
//
//PGDNMF::GradScatterMatrix PGDNMF::gradient_scatter(const Float2dArray& X, const Float2dArray& Z, const int& K, const Int1dArray& cls) {
//    int F = Z.rows();
//    int M = Z.cols();
//    int T = X.cols();
//    Float1dArray meanX = col_mean(X);
//    Float1dArray meanZ = col_mean(Z);
//    Float2dArray meanXcls = zeros(K, T);
//    Float2dArray meanZcls = zeros(K, M);
//    Int1dArray cnt(K);
//    for (int i = 0; i < F; ++i) {
//        int c = cls(i);
//        meanXcls(c, ALL) += X(i, ALL);
//        meanZcls(c, ALL) += Z(i, ALL);
//        cnt(c) += 1;
//    }
//    for (int c = 0; c < K; ++c) {
//        meanXcls(c, ALL) /= cnt(c);
//        meanZcls(c, ALL) /= cnt(c);
//    }
//    PGDNMF::GradScatterMatrix r;
//    r.grad_sw.resize(F, M);
//    r.gradd_sw.resize(M, M);
//    r.grad_sb.resize(F, M);
//    r.gradd_sb.resize(M, M);
//    r.grad_sw = 0;
//    //for (int i = 0; i < F; ++i)
//    //    for (int k = 0; k < M; ++k)
//    //        r.grad_sw(i, k) += x(i, 
//
//
//
//    r.grad_sb = 0;
//    r.gradd_sw = 0;
//    r.gradd_sb = 0;
//    return r;
//}
//
//PGDNMF::SubSolution PGDNMF::solve_Z(const Float2dArray& X, const Float2dArray& Z, const Float2dArray& H, const Int1dArray& cls, const double& tol, const size_t& maxiter) {
//}
//
//PGDNMF::SubSolution PGDNMF::solve_H(const Float2dArray& X, const Float2dArray& Z, const Float2dArray& H, const double& tol, const size_t& maxiter) {
//}
