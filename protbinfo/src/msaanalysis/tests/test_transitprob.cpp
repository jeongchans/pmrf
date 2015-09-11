#include <gtest/gtest.h>

#include "transitprob.h"

class TransitProbEstimatorTest : public testing::Test {
  public:
    TransitProbEstimatorTest() {
        uniform_rng.seed((unsigned int)time(0));
    }

    virtual void SetUp() {
        m_freq.resize(NUM_STATE_TYPE);
        d_freq.resize(NUM_STATE_TYPE);
        i_freq.resize(NUM_STATE_TYPE);
        for (size_t i = 0; i < NUM_STATE_TYPE; ++i) {
            m_freq(i) = uniform_rng.random();
            d_freq(i) = uniform_rng.random();
            i_freq(i) = uniform_rng.random();
        }
        d_freq(INSERT) = 0;
        i_freq(DELETE) = 0;
    }

    Float1dArray m_freq;
    Float1dArray d_freq;
    Float1dArray i_freq;
    ranlib::Uniform<double> uniform_rng;
};

TEST_F(TransitProbEstimatorTest, test_emd_transit_prob_estimator) {
    EMDTransitProbEstimator estimator;

    Float1dArray prob = estimator.estimate(m_freq, MATCH);
    ASSERT_EQ(NUM_STATE_TYPE, (size_t)prob.size());
    EXPECT_TRUE(allclose(1., blitz::sum(prob)));

    prob = estimator.estimate(d_freq, DELETE);
    ASSERT_EQ(NUM_STATE_TYPE, (size_t)prob.size());
    EXPECT_TRUE(allclose(1., blitz::sum(prob)));
    EXPECT_EQ(0, prob(INSERT));

    prob = estimator.estimate(i_freq, INSERT);
    ASSERT_EQ(NUM_STATE_TYPE, (size_t)prob.size());
    EXPECT_TRUE(allclose(1., blitz::sum(prob)));
    EXPECT_EQ(0, prob(DELETE));
}
