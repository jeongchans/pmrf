#include <gtest/gtest.h>

#include "seqweight.h"

class SeqWeightEstimatorTest : public testing::Test {
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

TEST_F(SeqWeightEstimatorTest, test_pb_seq_weight_estimator) {
    PBSeqWeightEstimator estimator;
    Float1dArray wt(4);
    wt = 0.267, 0.267, 0.200, 0.267;
    Float1dArray v = estimator.estimate(msa1);
    ASSERT_EQ(wt.size(), v.size());
    EXPECT_TRUE(allclose(wt, v)) << v;

    wt = 0.284, 0.230, 0.223, 0.263;
    v = estimator.estimate(msa2);
    EXPECT_TRUE(allclose(wt, v)) << v;
}
