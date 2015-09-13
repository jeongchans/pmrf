#include <gtest/gtest.h>

#include "mrfbuildcommand.h"

class MRFBuildCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFBuildCommandLine_Test, test_parse_param) {
    int argc = 5;
    char* argv[5] = {"pmrf", "build",
                     "aaa.afa", "-o", "aaa.mrf"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(AFASTA, cmd_line.opt.build_opt.msa_fmt);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_regul_param) {
    int argc = 15;
    char* argv[15] = {"pmrf", "build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--regnode-l2", "0",
                      "--l2-lambda1", "2.0",
                      "--regedge-l2", "1",
                      "--l2-lambda2", "3.0",
                      "--l2-scale", "0"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(false, cmd_line.opt.build_opt.parameterizer_opt.node_l2_regul);
    EXPECT_EQ(2.0, cmd_line.opt.build_opt.parameterizer_opt.node_l2_opt.lambda);

    EXPECT_EQ(true, cmd_line.opt.build_opt.parameterizer_opt.edge_l2_regul);
    EXPECT_EQ(3.0, cmd_line.opt.build_opt.parameterizer_opt.edge_l2_opt.lambda);
    EXPECT_EQ(false, cmd_line.opt.build_opt.parameterizer_opt.edge_l2_opt.sc);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_input_param) {
    int argc = 9;
    char* argv[9] = {"pmrf", "build",
                     "aaa.afa", "-o", "aaa.mrf",
                     "--msa", "a3m",
                     "--edge", "edges.txt"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(A3M, cmd_line.opt.build_opt.msa_fmt);

    EXPECT_EQ("edges.txt", cmd_line.opt.build_opt.eidx_filename);
}
