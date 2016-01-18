#include <gtest/gtest.h>

#include "mrfmaincommand.h"

class MRFMainCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFMainCommandLine_Test, test_invalid_subcmd) {
    int argc = 2;
    char* argv[2] = {"pmrf", "undetermined"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_FALSE(cmd_line.is_valid());
}

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_help) {
    int argc = 2;
    char* argv[2] = {"pmrf", "help"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(HELP, cmd_line.subcmd);
}

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_build) {
    int argc = 5;
    char* argv[5] = {"pmrf", "build",
                     "aaa.afa", "-o", "aaa.mrf"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(BUILD, cmd_line.subcmd);
}

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_infer) {
    int argc = 4;
    char* argv[4] = {"pmrf", "infer",
                     "aaa.mrf", "aaa.fa"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(INFER, cmd_line.subcmd);
}
