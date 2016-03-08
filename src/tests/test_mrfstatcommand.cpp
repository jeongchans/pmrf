#include <gtest/gtest.h>

#include "mrfstatcommand.h"

class MRFStatCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFStatCommandLine_Test, test_parse_param) {
    int argc = 3;
    char* argv[3] = {"pmrf", "stat",
                     "aaa.mrf"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.mrf", cmd_line.opt.mrf_filename);

    EXPECT_EQ(STATMODE_PAIR, cmd_line.opt.stat_opt.mode);
    EXPECT_EQ(STATCORR_APC, cmd_line.opt.stat_opt.corr);
}

TEST_F(MRFStatCommandLine_Test, test_parse_opt_param) {
    int argc = 7;
    char* argv[7] = {"pmrf", "stat",
                     "aaa.mrf",
                     "--mode", "pos",
                     "--corr", "2"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.mrf", cmd_line.opt.mrf_filename);

    EXPECT_EQ(STATMODE_POS, cmd_line.opt.stat_opt.mode);
    EXPECT_EQ(STATCORR_NCPS, cmd_line.opt.stat_opt.corr);
}
