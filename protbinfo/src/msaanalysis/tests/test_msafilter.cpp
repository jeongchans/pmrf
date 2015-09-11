#include <gtest/gtest.h>

#include "msafilter.h"

class MSAFilterTest : public testing::Test {
  public:
    virtual void SetUp() {
        msa.push_back("--VGS");
        msa.push_back("-FDGF");
        msa.push_back("GYDG-");
        msa.push_back("GYQG-");
    }

    AminoAcid abc;
    vector<string> msa;
};

TEST_F(MSAFilterTest, test_terminal_gap_remover) {
    TerminalGapRemover filter(abc, 0.3);
    vector<string> filt = filter.filter(msa);
    ASSERT_EQ(msa.size(), filt.size());
    EXPECT_EQ("-VG", filt[0]);
    EXPECT_EQ("FDG", filt[1]);
    EXPECT_EQ("YDG", filt[2]);
    EXPECT_EQ("YQG", filt[3]);
}
