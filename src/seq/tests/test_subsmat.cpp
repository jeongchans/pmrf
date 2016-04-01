#include <gtest/gtest.h>

#include "subsmat.h"

#include <random/uniform.h>

using ranlib::Uniform;

class SubstitutionMatrixTest : public testing::Test {
  public:
    AminoAcid abc;
};

TEST_F(SubstitutionMatrixTest, test_substitution_matrix) {
    std::string s = abc.get_canonical();
    std::map<SymbolPair, double> score;
    MatrixXf m(s.size(), s.size());
    for (int i = 0; i < (int)s.size(); ++i) {
        for (int j = i; j < (int)s.size(); ++j) {
            double val = i + j;
            score.insert(make_pair(SymbolPair(s[i], s[j]), val));
            m(i, j) = m(j, i) = val;
        }
    }
    SubstitutionMatrix subs_mat(score);
    EXPECT_TRUE(m.matrix() == subs_mat.get_array(abc).matrix());
    for (int i = 0; i < (int)s.size(); ++i) {
        for (int j = i; j < (int)s.size(); ++j) {
            double val = i + j;
            EXPECT_EQ(val, subs_mat.get_value(s[i], s[j]));
            EXPECT_EQ(val, subs_mat.get_value(s[j], s[i]));
        }
    }
}

TEST_F(SubstitutionMatrixTest, test_blosum62) {
    BLOSUM62Matrix blosum62;
    MatrixXf m(20, 20);
    //    A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
    m <<  4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2,  // A
          0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,  // C
         -2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3,  // D
         -1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2,  // E
         -2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3,  // F
          0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3,  // G
         -2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2,  // H
         -1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1,  // I
         -1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2,  // K
         -1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1,  // L
         -1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1,  // M
         -2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2,  // N
         -1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3,  // P
         -1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1,  // Q
         -1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2,  // R
          1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2,  // S
          0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2,  // T
          0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1,  // V
         -3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2,  // W
         -2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7;  // Y
    EXPECT_TRUE(m.matrix() == blosum62.get_array(abc).matrix());
    EXPECT_EQ(4, blosum62.get_value('A', 'A'));
    EXPECT_EQ(-3, blosum62.get_value('C', 'D'));
    EXPECT_EQ(-3, blosum62.get_value('D', 'C'));
}

class TargetProbEstimatorTest : public testing::Test {
  public:
    TargetProbEstimatorTest() {
        rng.seed((unsigned int)time(0));
    }

    virtual void SetUp() {
        lamb = 0.3;
        prob.resize(20, 20);
        for (int i = 0; i < 20; ++i) {
            prob(i, i) = rng.random() + 0.5;
            for (int j = i + 1; j < 20; ++j)
                prob(i, j) = prob(j, i) = rng.random();
        }
        prob /= prob.sum();
        bgfreq.resize(20);
        for (int i = 0; i < 20; ++i) bgfreq(i) = rng.random();
        bgfreq /= bgfreq.sum();
        subsmat.resize(prob.rows(), prob.cols());
        subsmat = 1. / lamb * prob.cwiseQuotient(bgfreq * bgfreq.transpose()).unaryExpr(&log);
    }

    Uniform<double> rng;
    VectorXf bgfreq;
    double lamb;
    MatrixXf prob;
    MatrixXf subsmat;
};

TEST_F(TargetProbEstimatorTest, test_target_prob_estimator_given_bg) {
    TargetProbEstimatorGivenBG estimator(bgfreq);
    pair<double, MatrixXf> p = estimator.probify(subsmat);
    EXPECT_TRUE(allclose(lamb, p.first)) << p.first;
    ASSERT_EQ(prob.rows(), p.second.rows());
    ASSERT_EQ(prob.cols(), p.second.cols());
    EXPECT_TRUE(allclose(prob, p.second));
}
