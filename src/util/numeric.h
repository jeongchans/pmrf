#ifndef _PROTBINFO_NUMERIC_H_
#define _PROTBINFO_NUMERIC_H_

#ifndef FLOAT_PRECISION
#define FLOAT_PRECISION     32
#endif

#include <Eigen/Core>

//#include <cmath>
//#include <ctime>
//#include <iomanip>

//#include <blitz/array.h>
//#include <random/normal.h>

#if     FLOAT_PRECISION == 32
typedef float FloatType;
#elif   FLOAT_PRECISION == 64
typedef double FloatType;
#endif

using Eigen::Dynamic;

//using blitz::abs;
//using blitz::all;
//using blitz::Range;

//static Range ALL = Range::all();

typedef Eigen::Array<FloatType, Dynamic, 1> Float1dArray;
typedef Eigen::Array<FloatType, Dynamic, Dynamic> Float2dArray;

//typedef blitz::Array<int, 1> Int1dArray;
//typedef blitz::Array<FloatType, 1> Float1dArray;
//typedef blitz::Array<FloatType, 2> Float2dArray;
//typedef blitz::Array<FloatType, 3> Float3dArray;
//typedef blitz::Array<FloatType, 4> Float4dArray;

inline bool allclose(const Float2dArray& x, const Float2dArray& y, const FloatType& atol=1e-3) {
    return x.isApprox(y, atol);
}

inline bool allclose(const Float1dArray& x, const Float1dArray& y, const FloatType& atol=1e-3) {
    return x.isApprox(y, atol);
}

inline bool allclose(const FloatType& x, const FloatType& y, const FloatType& atol=1e-3) {
    return fabs(x - y) < atol;
}

inline FloatType sum(const Float2dArray& x) {
    return x.sum();
}

inline FloatType mean(const Float1dArray& x) {
    return x.mean();
}

//inline FloatType norm(const Float1dArray& x) {
//    return sqrt(blitz::sum(blitz::pow2(x)));
//}

inline FloatType norm(const Float2dArray& x) {
    return x.matrix().norm();
}

//inline Float2dArray dot(const Float2dArray& x, const Float2dArray& y) {
//    Float2dArray r(x.rows(), y.cols());
//    blitz::firstIndex i;
//    blitz::secondIndex j;
//    blitz::thirdIndex k;
//    r = blitz::sum(x(i, k) * y(k, j), k);
//    return r;
//}

//inline Float1dArray dot(const Float2dArray& x, const Float1dArray& y) {
//    Float1dArray r(x.rows());
//    blitz::firstIndex i;
//    blitz::secondIndex j;
//    r = blitz::sum(x(i, j) * y(j), j);
//    return r;
//}

inline FloatType max(const Float2dArray& x) {
    return x.maxCoeff();
}

inline Float1dArray row_max(const Float2dArray& x) {
    return x.rowwise().maxCoeff();
}

//inline Float1dArray row_sum(const Float2dArray& x) {
//    Float1dArray r(x.rows());
//    blitz::secondIndex j;
//    r = blitz::sum(x, j);
//    return r;
//}

//inline Float1dArray col_sum(const Float2dArray& x) {
//    Float1dArray r(x.cols());
//    blitz::firstIndex i;
//    blitz::secondIndex j;
//    r = blitz::sum(x(j, i), j);
//    return r;
//}

//inline Float1dArray row_mean(const Float2dArray& x) {
//    Float1dArray r(x.rows());
//    blitz::secondIndex j;
//    r = blitz::mean(x, j);
//    return r;
//}

//inline Float1dArray col_mean(const Float2dArray& x) {
//    Float1dArray r(x.cols());
//    blitz::firstIndex i;
//    blitz::secondIndex j;
//    r = blitz::mean(x(j, i), j);
//    return r;
//}

//inline Float1dArray row_stdev(const Float2dArray& x) {
//    using blitz::pow2;
//    Float2dArray x2(x.shape());
//    x2 = pow2(x);
//    Float1dArray mn = row_mean(x);
//    Float1dArray mn2 = row_mean(x2);
//    Float1dArray r(x.rows());
//    r = sqrt(mn2 - pow2(mn));
//    return r;
//}

//inline Float1dArray col_stdev(const Float2dArray& x) {
//    using blitz::pow2;
//    Float2dArray x2(x.shape());
//    x2 = pow2(x);
//    Float1dArray mn = col_mean(x);
//    Float1dArray mn2 = col_mean(x2);
//    Float1dArray r(x.cols());
//    r = sqrt(mn2 - pow2(mn));
//    return r;
//}

//inline FloatType stdev(const Float2dArray& x) {
//    using blitz::mean;
//    using blitz::pow2;
//    return sqrt(mean(pow2(x)) - pow2(mean(x)));
//}

inline Float2dArray outer(const Float1dArray& x, const Float1dArray& y) {
    return (x.matrix() * y.matrix().transpose()).array();
    //return x * y.transpose();
}

inline Float2dArray transpose(const Float2dArray& x) {
    return x.transpose();
}

//inline Float2dArray concatenate(const Float2dArray& x, const Float2dArray& y) {
//    if (x.cols() != y.cols()) throw;
//    int cols = x.cols();
//    int n1 = x.rows();
//    int n2 = y.rows();
//    Float2dArray r(n1 + n2, cols);
//    r(Range(blitz::fromStart, n1 - 1), ALL) = x.copy();
//    r(Range(n1, blitz::toEnd), ALL) = y.copy();
//    return r;
//}

//inline Float1dArray flatten(const Float2dArray& x) {
//    int rows = x.rows();
//    int cols = x.cols();
//    Float1dArray r(rows * cols);
//    for (int i = 0; i < rows; ++i)
//        r(Range(i * cols, (i + 1) * cols - 1)) = x(i, ALL).copy();
//    return r;
//}

inline void scale(Float1dArray& v, const FloatType& sc=1.0) {
    FloatType s = v.sum();
    if (s != 0) {
        v /= s;
        v *= sc;
    }
}

//inline void scale(Float2dArray& v, const FloatType& sc=1.0) {
//    FloatType s = blitz::sum(v);
//    if (s != 0) {
//        v /= s;
//        v *= sc;
//    }
//}

inline Float2dArray pow2(const Float2dArray& m) {
    return m.square();
}

inline Float1dArray pow2(const Float1dArray& v) {
    return v.square();
}

//inline Float1dArray exp(const Float1dArray& v) {
//    int n = v.size();
//    Float1dArray r(n);
//    for (int i = 0; i < n; ++i) r(i) = exp(v(i));
//    return r;
//}

//inline Float2dArray exp(const Float2dArray& m) {
//    int rows = m.rows();
//    int cols = m.cols();
//    Float2dArray r(rows, cols);
//    for (int i = 0; i < rows; ++i)
//        for (int j = 0; j < cols; ++j)
//            r(i, j) = exp(m(i, j));
//    return r;
//}

//inline Float1dArray log(const Float1dArray& v) {
//    int n = v.size();
//    Float1dArray r(n);
//    for (int i = 0; i < n; ++i) r(i) = log(v(i));
//    return r;
//}

//inline Float2dArray log(const Float2dArray& m) {
//    int rows = m.rows();
//    int cols = m.cols();
//    Float2dArray r(rows, cols);
//    for (int i = 0; i < rows; ++i)
//        for (int j = 0; j < cols; ++j)
//            r(i, j) = log(m(i, j));
//    return r;
//}

//inline Float1dArray log2(const Float1dArray& v) {
//    int n = v.size();
//    Float1dArray r(n);
//    for (int i = 0; i < n; ++i) r(i) = log2(v(i));
//    return r;
//}

//inline Float2dArray log2(const Float2dArray& m) {
//    int rows = m.rows();
//    int cols = m.cols();
//    Float2dArray r(rows, cols);
//    for (int i = 0; i < rows; ++i)
//        for (int j = 0; j < cols; ++j)
//            r(i, j) = log2(m(i, j));
//    return r;
//}

//inline FloatType stddev(const Float1dArray& v) {
//    FloatType mn = blitz::mean(v);
//    return sqrt(blitz::mean(blitz::pow2(v - mn)));
//}

inline Float1dArray zeros(const size_t& n) {
    return Float1dArray::Zero(n);
}

inline Float2dArray zeros(const size_t& n, const size_t& m) {
    return Float2dArray::Zero(n, m);
}

inline void init_rng() {
    static bool seeded = false;
    if (!seeded) {
        srand((unsigned int) time(0));
        seeded = true;
    }
}

inline Float2dArray randu(const size_t& n, const size_t& m) {
    init_rng();
    return (Float2dArray::Random(n, m) + 1.) * 0.5;
}

inline Float1dArray randu(const size_t& n) {
    init_rng();
    return (Float1dArray::Random(n) + 1.) * 0.5;
}

inline FloatType randu() {
    return randu(1)(0);
}

const FloatType TWO_PI = 2. * 3.14159265358979323846;

inline Float2dArray randn(const size_t& n, const size_t& m) {
    Float2dArray u1 = randu(n, m);
    Float2dArray u2 = randu(n, m);
    return sqrt(-2. * log(u1)) * cos(TWO_PI * u2);
}

inline Float1dArray randn(const size_t& n) {
    Float1dArray u1 = randu(n);
    Float1dArray u2 = randu(n);
    return sqrt(-2. * log(u1)) * cos(TWO_PI * u2);
}

inline FloatType randn() {
    return randn(1)(0);
}

//static inline void export_array_elem(const FloatType& x, std::ostream& os) {
//    os.unsetf(std::ios::floatfield);
//    os << std::setprecision(6) << x;
//}

//inline void export_array(const Float1dArray& x, std::ostream& os) {
//    using std::endl;
//    int n = x.size();
//    os << n << endl;
//    for (int i = 0; i < n; ++i) {
//        export_array_elem(x(i), os);
//        if (i != n - 1) os << "\t";
//    }
//    os << endl;
//}

//inline void export_array(const Float2dArray& x, std::ostream& os) {
//    using std::endl;
//    int rows = x.rows();
//    int cols = x.cols();
//    os << rows << " x " << cols << endl;
//    for (int i = 0; i < rows; ++i) {
//        for (int j = 0; j < cols; ++j) {
//            export_array_elem(x(i, j), os);
//            if (j != cols - 1) os << "\t";
//        }
//        os << endl;
//    }
//}

//inline void import_array(Float1dArray& x, std::istream& is) {
//    int n;
//    is >> n;
//    x.resize(n);
//    for (int i = 0; i < n; ++i)
//        is >> x(i);
//}

//inline void import_array(Float2dArray& x, std::istream& is) {
//    std::string dummy;
//    int rows, cols;
//    is >> rows >> dummy >> cols;
//    x.resize(rows, cols);
//    for (int i = 0; i < rows; ++i)
//        for (int j = 0; j < cols; ++j)
//            is >> x(i, j);
//}

//inline std::vector<const FloatType*> raw_data(const Float2dArray& x) {
//    const int rows = x.rows();
//    const int cols = x.cols();
//    const FloatType* ptr = x.data();
//    std::vector<const FloatType*> r(rows);
//    for (int i = 0; i < rows; ++i)
//        r[i] = &ptr[i * cols];
//    return r;
//}

inline FloatType calc_kld(const Float1dArray f, const Float1dArray p) {
    FloatType d = 0.;
    for (size_t t = 0; t < f.size(); ++t)
        if (f(t) > 0) d += f(t) * log(f(t) / p(t));
    return d;
}

#endif
