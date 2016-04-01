#include <gtest/gtest.h>

#include "probdistrib.h"

#include "bgfreq.h"
#include "subsmat.h"

class EmitProbEstimatorTest : public testing::Test {
  public:
    virtual void SetUp() {
        AminoAcid abc;
        bgfreq = RobinsonBgFreq().get_array(abc);
        scoremat = BLOSUM62Matrix().get_array(abc);
        freq = VectorXf::Random(20).cwiseAbs();
    }

    VectorXf bgfreq;
    MatrixXf scoremat;
    VectorXf freq;
};

TEST_F(EmitProbEstimatorTest, test_smm_emit_prob_estimator) {
    SMMEmitProbEstimator estimator(bgfreq, scoremat);
    VectorXf prob = estimator.estimate(freq);
    ASSERT_EQ(freq.size(), prob.size());
    EXPECT_FLOAT_EQ(1., prob.sum());
}
