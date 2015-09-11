#include <gtest/gtest.h>

#include "emitprob.h"

#include "seq/bgfreq.h"
#include "seq/subsmat.h"

class EmitProbEstimatorTest : public testing::Test {
  public:
    EmitProbEstimatorTest() {
        uniform_rng.seed((unsigned int)time(0));
    }

    virtual void SetUp() {
        AminoAcid abc;
        bgfreq.resize(20);
        bgfreq = RobinsonBgFreq().get_array(abc);
        scoremat.resize(20, 20);
        scoremat = BLOSUM62Matrix().get_array(abc);
        freq.resize(20);
        for (int i = 0; i < 20; ++i) freq(i) = uniform_rng.random();
    }

    Float1dArray bgfreq;
    Float2dArray scoremat;
    Float1dArray freq;
    ranlib::Uniform<double> uniform_rng;
};

TEST_F(EmitProbEstimatorTest, test_smm_emit_prob_estimator) {
    SMMEmitProbEstimator estimator(bgfreq, scoremat);
    Float1dArray prob = estimator.estimate(freq);
    ASSERT_EQ(freq.size(), prob.size());
    EXPECT_TRUE(allclose(1., blitz::sum(prob)));
}
