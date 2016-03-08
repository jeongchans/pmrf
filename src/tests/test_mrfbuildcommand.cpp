#include <gtest/gtest.h>

#include "mrfbuildcommand.h"

class MRFBuildCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFBuildCommandLine_Test, test_parse_param_default) {
    int argc = 5;
    char* argv[5] = {"pmrf", "build",
                     "aaa.afa", "-o", "aaa.mrf"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(AFASTA, cmd_line.opt.build_opt.msa_fmt);

    EXPECT_EQ(POSITION_BASED, cmd_line.opt.build_opt.msa_analyzer_opt.seq_wt);
    EXPECT_EQ(NO_EFFNUM, cmd_line.opt.build_opt.msa_analyzer_opt.eff_num);

    EXPECT_EQ(RegulMethod::RegulMethod::L2, cmd_line.opt.build_opt.parameterizer_opt.regul);
    EXPECT_EQ(0.01, cmd_line.opt.build_opt.parameterizer_opt.l2_opt.lambda1);
    EXPECT_EQ(0.2, cmd_line.opt.build_opt.parameterizer_opt.l2_opt.lambda2);
    EXPECT_EQ(true, cmd_line.opt.build_opt.parameterizer_opt.l2_opt.sc);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_l2_regul_param) {
    int argc = 13;
    char* argv[13] = {"pmrf", "build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--regul", "1",
                      "--regnode-lambda", "15.0",
                      "--regedge-lambda", "3.0",
                      "--regedge-scale", "0"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(RegulMethod::RegulMethod::L2, cmd_line.opt.build_opt.parameterizer_opt.regul);
    EXPECT_EQ(15.0, cmd_line.opt.build_opt.parameterizer_opt.l2_opt.lambda1);
    EXPECT_EQ(3.0, cmd_line.opt.build_opt.parameterizer_opt.l2_opt.lambda2);
    EXPECT_EQ(false, cmd_line.opt.build_opt.parameterizer_opt.l2_opt.sc);
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

TEST_F(MRFBuildCommandLine_Test, test_parse_preproc_param) {
    int argc = 9;
    char* argv[9] = {"pmrf", "build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--seqwt", "0",
                      "--effnum", "1"};
    MRFBuildCommandLine cmd_line(argc, argv);

    EXPECT_EQ(NO_WEIGHT, cmd_line.opt.build_opt.msa_analyzer_opt.seq_wt);
    EXPECT_EQ(EXP_ENTROPY, cmd_line.opt.build_opt.msa_analyzer_opt.eff_num);
}
