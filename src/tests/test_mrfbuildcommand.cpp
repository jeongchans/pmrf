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

TEST_F(MRFBuildCommandLine_Test, test_parse_node_regul_param) {
    int argc = 11;
    char* argv[11] = {"pmrf", "build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--regnode", "2",
                      "--regnode-lambda", "15.0",
                      "--regedge", "0"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(NodeRegulMethod::PSSM, cmd_line.opt.build_opt.parameterizer_opt.node_regul);
    EXPECT_EQ(15.0, cmd_line.opt.build_opt.parameterizer_opt.node_pssm_opt.lambda);

    EXPECT_EQ(EdgeRegulMethod::NONE, cmd_line.opt.build_opt.parameterizer_opt.edge_regul);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_edge_regul_param) {
    int argc = 9;
    char* argv[9] = {"pmrf", "build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--regedge-lambda", "3.0",
                      "--regedge-scale", "0"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(EdgeRegulMethod::L2, cmd_line.opt.build_opt.parameterizer_opt.edge_regul);
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
