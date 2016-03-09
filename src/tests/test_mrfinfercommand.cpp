#include <gtest/gtest.h>

#include "mrfinfercommand.h"

class MRFInferCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFInferCommandLine_Test, test_parse_param) {
    int argc = 3;
    char* argv[3] = {"infer",
                     "aaa.mrf", "aaa.fa"};
    MRFInferCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.mrf", cmd_line.opt.mrf_filename);
    EXPECT_EQ("aaa.fa", cmd_line.opt.seq_filename);
}
