#ifndef _PROTBINFO_NUMERIC_H_
#define _PROTBINFO_NUMERIC_H_

#include <random>

#include <eigen3/Eigen/Core>

typedef float FloatType;

//using Eigen::Dynamic;
using Eigen::VectorXf;

using Eigen::MatrixXi;
using Eigen::MatrixXf;
using Eigen::MatrixXd;

//typedef Eigen::Array<FloatType, Dynamic, 1> Float1dArray;
//typedef Eigen::Array<FloatType, Dynamic, Dynamic> Float2dArray;

//inline bool allclose(const Float2dArray& x, const Float2dArray& y, const FloatType& atol=1e-3) {
//    return x.isApprox(y, atol);
//}

//inline bool allclose(const Float1dArray& x, const Float1dArray& y, const FloatType& atol=1e-3) {
//    return x.isApprox(y, atol);
//}

inline bool allclose(const FloatType& x, const FloatType& y, const FloatType& atol=1e-3) {
    return fabs(x - y) < atol;
}

inline bool allclose(const VectorXf& x, const VectorXf& y, const float& atol=1e-3) {
    return x.isApprox(y, atol);
}

inline bool allclose(const MatrixXf& x, const MatrixXf& y, const float& atol=1e-3) {
    return x.isApprox(y, atol);
}

inline float square(const float x) {
    return x * x;
}

//inline FloatType sum(const Float2dArray& x) {
//    return x.sum();
//}

//inline FloatType mean(const Float1dArray& x) {
//    return x.mean();
//}

//inline FloatType norm(const Float2dArray& x) {
//    return x.matrix().norm();
//}

//inline FloatType max(const Float2dArray& x) {
//    return x.maxCoeff();
//}

//inline Float1dArray row_max(const Float2dArray& x) {
//    return x.rowwise().maxCoeff();
//}

//inline Float2dArray outer(const Float1dArray& x, const Float1dArray& y) {
//    return (x.matrix() * y.matrix().transpose()).array();
//    //return x * y.transpose();
//}

//inline Float2dArray transpose(const Float2dArray& x) {
//    return x.transpose();
//}

//inline void scale(Float1dArray& v, const FloatType& sc=1.0) {
//    FloatType s = v.sum();
//    if (s != 0) {
//        v /= s;
//        v *= sc;
//    }
//}

//inline Float2dArray pow2(const Float2dArray& m) {
//    return m.square();
//}

//inline Float1dArray pow2(const Float1dArray& v) {
//    return v.square();
//}

//inline Float1dArray zeros(const size_t& n) {
//    return Float1dArray::Zero(n);
//}

//inline Float2dArray zeros(const size_t& n, const size_t& m) {
//    return Float2dArray::Zero(n, m);
//}

inline void init_rng() {
    static bool seeded = false;
    if (!seeded) {
        srand((unsigned int) time(0));
        seeded = true;
    }
}

//inline Float2dArray randu(const size_t& n, const size_t& m) {
//    init_rng();
//    return (Float2dArray::Random(n, m) + 1.) * 0.5;
//}

//inline Float1dArray randu(const size_t& n) {
//    init_rng();
//    return (Float1dArray::Random(n) + 1.) * 0.5;
//}

//inline FloatType randu() {
//    return randu(1)(0);
//}

//const float TWO_PI = 2. * 3.14159265358979323846;

inline float randn() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0., 1.);
    return d(gen);
}

inline VectorXf randn_vector(const size_t n) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0., 1.);
    VectorXf v(n);
    for (size_t i = 0; i < n; ++i) v(i) = d(gen);
    return v;
}

inline MatrixXf randn_matrix(const size_t n, const size_t m) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0., 1.);
    MatrixXf x(n, m);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < m; ++j) x(i, j) = d(gen);
    return x;
}

//inline Float2dArray randn(const size_t& n, const size_t& m) {
//    Float2dArray u1 = randu(n, m);
//    Float2dArray u2 = randu(n, m);
//    return sqrt(-2. * log(u1)) * cos(TWO_PI * u2);
//}

//inline Float1dArray randn(const size_t& n) {
//    Float1dArray u1 = randu(n);
//    Float1dArray u2 = randu(n);
//    return sqrt(-2. * log(u1)) * cos(TWO_PI * u2);
//}

//inline FloatType randn() {
//    return randn(1)(0);
//}

inline float kld_elem(const float p, const float q) {
    if (p > 0) return p * log(p / q);
    return 0.;
}

inline float kld(const VectorXf p, const VectorXf q) {
    return p.binaryExpr(q, &kld_elem).sum();
}

#endif
