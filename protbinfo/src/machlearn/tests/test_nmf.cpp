#include <gtest/gtest.h>

#include "nmf.h"

class NMFTest : public testing::Test {
  protected:
    virtual void SetUp() {
        n = 100;
        m = 3;
        r = 2;
        W0.resize(m, r);
        H0.resize(n, r);
        W0 = blitz::abs(randn(m, r));
        H0 = blitz::abs(randn(n, r));
        V.resize(n, m);
        V = transpose(dot(W0, transpose(H0)));
    }

    int n;  // number of samples
    int m;  // number of original fatures
    int r;  // number of basis
    Float2dArray W0;
    Float2dArray H0;
    Float2dArray V;
};

TEST_F(NMFTest, test_nonnegative_input) {
    EXPECT_TRUE(all(V >= 0));
    EXPECT_EQ(n, V.rows());
    EXPECT_EQ(m, V.cols());
}

TEST_F(NMFTest, test_nmf) {
    Float2dArray W, H;
    NMF nmf;

    nmf.factorize(V, r, W, H);
    ASSERT_TRUE(W0.rows() == W.rows() && W0.cols() == W.cols()) << W.shape();
    EXPECT_TRUE(all(W >= 0));
    ASSERT_TRUE(H0.rows() == H.rows() && H0.cols() == H.cols()) << H.shape();
    EXPECT_TRUE(all(H >= 0));
    EXPECT_TRUE(allclose(V, transpose(dot(W, transpose(H)))));

    EXPECT_TRUE(allclose(H, nmf.transform(V, W)));

    Float1dArray x = V(0, ALL);
    Float1dArray y = nmf.transform(x, W);
    ASSERT_EQ(r, y.size());
    EXPECT_TRUE(allclose(H(0, ALL), y));
}

//class PGDNMFTest : public testing::Test {
//  protected:
//    virtual void SetUp() {
//        using blitz::fromStart;
//        using blitz::toEnd;
//        n = 100;
//        m = 3;
//        r = 2;
//        int num_cls1 = 50;
//        W0.resize(m, r);
//        H0.resize(n, r);
//        W0 = blitz::abs(randn(m, r));
//        H0 = blitz::abs(randn(n, r));
//        H0(Range(fromStart, num_cls1 - 1)) += 3.;
//        V.resize(n, m);
//        V = transpose(dot(W0, transpose(H0)));
//        cls.resize(n);
//        cls(Range(fromStart, num_cls1- - 1)) = 0;
//        cls(Range(num_cls1, toEnd)) = 1;
//        cls_wt.resize(2);
//    }
//
//    int n;  // number of samples
//    int m;  // number of original fatures
//    int r;  // number of basis
//    Float2dArray W0;   // basis matrix
//    Float2dArray H0;   // weight matrix
//    Float2dArray V;
//    Int1dArray cls;  // class definition
//    Float1dArray cls_wt;   // class weight
//};
//
//TEST_F(PGDNMFTest, test_nonnegative_input) {
//    EXPECT_TRUE(all(V >= 0));
//    EXPECT_EQ(n, V.rows());
//    EXPECT_EQ(m, V.cols());
//}
//
//TEST_F(PGDNMFTest, test_factorize) {
//    Float2dArray W, H;
//    PGDNMF nmf;
//
//    NMFMatrix m = nmf.factorize(V, r, 2, cls, cls_wt);
//    W = m.basis;
//    H = m.weight;
//    ASSERT_EQ(W0.rows(), W.rows());
//    ASSERT_EQ(W0.cols(), W.cols());
//    EXPECT_TRUE(all(W >= 0));
//    ASSERT_EQ(H0.rows(), H.rows());
//    ASSERT_EQ(H0.cols(), H.cols());
//    EXPECT_TRUE(all(H >= 0));
//
//    EXPECT_TRUE(allclose(V, transpose(dot(W, transpose(H)))));
//
//    EXPECT_TRUE(allclose(V, nmf.transform(V, W)));
//
//    Float1dArray x = V(0, ALL);
//    Float1dArray y = nmf.transform(x, W);
//    ASSERT_EQ(r, y.size());
//    EXPECT_TRUE(allclose(H(0, ALL), y));
//}
