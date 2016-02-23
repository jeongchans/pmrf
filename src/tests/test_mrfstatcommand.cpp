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
}
