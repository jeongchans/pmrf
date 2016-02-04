#include <gtest/gtest.h>

#include "effseqnum.h"

class EffSeqNumEstimatorTest : public testing::Test {
  public:
    virtual void SetUp() {
        msa1.push_back("GYVGS");
        msa1.push_back("GFDGF");
        msa1.push_back("GYDGF");
        msa1.push_back("GYQGG");

        msa2.push_back("--VGS");
        msa2.push_back("-FDGF");
        msa2.push_back("GYDGF");
        msa2.push_back("GYQG-");
    }

    AminoAcid abc;
    vector<string> msa1;    // MSA without gaps
    vector<string> msa2;    // MSA with gaps
};

TEST_F(EffSeqNumEstimatorTest, test_rt_eff_seq_num_estimator) {
    RTEffSeqNumEstimator estimator(abc);
    EXPECT_FLOAT_EQ(2.0, estimator.estimate(msa1));
    EXPECT_FLOAT_EQ(1.8, estimator.estimate(msa2));
}

TEST_F(EffSeqNumEstimatorTest, test_exp_entropy_eff_seq_num_estimator) {
    NullSeqWeightEstimator seq_weight_estimator;
    ExpEntropyEffSeqNumEstimator estimator(abc, &seq_weight_estimator);
    //EXPECT_TRUE(allclose(2.1431, estimator.estimate(msa1))) << estimator.estimate(msa1);  // base=2, range=[1, 75.33]
    EXPECT_TRUE(allclose(1.6961, estimator.estimate(msa1))) << estimator.estimate(msa1);    // base=e, range=[1, 20]
}