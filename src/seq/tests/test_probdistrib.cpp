#include <gtest/gtest.h>

#include "probdistrib.h"

#include "bgfreq.h"
#include "subsmat.h"

class EmitProbEstimatorTest : public testing::Test {
  public:
    virtual void SetUp() {
        AminoAcid abc;
        bgfreq.resize(20);
        bgfreq = RobinsonBgFreq().get_array(abc);
        scoremat.resize(20, 20);
        scoremat = BLOSUM62Matrix().get_array(abc);
        freq.resize(20);
        for (int i = 0; i < 20; ++i) freq(i) = randu();
    }

    Float1dArray bgfreq;
    Float2dArray scoremat;
    Float1dArray freq;
};

TEST_F(EmitProbEstimatorTest, test_smm_emit_prob_estimator) {
    SMMEmitProbEstimator estimator(bgfreq, scoremat);
    Float1dArray prob = estimator.estimate(freq);
    ASSERT_EQ(freq.size(), prob.size());
    EXPECT_FLOAT_EQ(1., sum(prob));
}
