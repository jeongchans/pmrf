#include <gtest/gtest.h>

#include "command.h"

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

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_stat) {
    int argc = 3;
    char* argv[3] = {"pmrf", "stat",
                     "aaa.mrf"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(STAT, cmd_line.subcmd);
}

class MRFBuildCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFBuildCommandLine_Test, test_parse_param_default) {
    int argc = 4;
    char* argv[4] = {"build",
                     "aaa.afa", "-o", "aaa.mrf"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(AFASTA, cmd_line.opt.msa_fmt);

    EXPECT_EQ(POSITION_BASED, cmd_line.opt.msa_analyzer_opt.seq_wt);
    EXPECT_EQ(NO_EFFNUM, cmd_line.opt.msa_analyzer_opt.eff_num);

    EXPECT_EQ(RegulMethod::RegulMethod::L2, cmd_line.opt.parameterizer_opt.regul);
    EXPECT_EQ(0.01, cmd_line.opt.parameterizer_opt.l2_opt.lambda1);
    EXPECT_EQ(0.2, cmd_line.opt.parameterizer_opt.l2_opt.lambda2);
    EXPECT_EQ(true, cmd_line.opt.parameterizer_opt.l2_opt.sc);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_l2_regul_param) {
    int argc = 12;
    char* argv[12] = {"build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--regul", "1",
                      "--regnode-lambda", "15.0",
                      "--regedge-lambda", "3.0",
                      "--regedge-scale", "0"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(RegulMethod::RegulMethod::L2, cmd_line.opt.parameterizer_opt.regul);
    EXPECT_EQ(15.0, cmd_line.opt.parameterizer_opt.l2_opt.lambda1);
    EXPECT_EQ(3.0, cmd_line.opt.parameterizer_opt.l2_opt.lambda2);
    EXPECT_EQ(false, cmd_line.opt.parameterizer_opt.l2_opt.sc);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_input_param) {
    int argc = 8;
    char* argv[8] = {"build",
                     "aaa.afa", "-o", "aaa.mrf",
                     "--msa", "a3m",
                     "--edge", "edges.txt"};
    MRFBuildCommandLine cmd_line(argc, argv);

    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("aaa.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(A3M, cmd_line.opt.msa_fmt);

    EXPECT_EQ("edges.txt", cmd_line.opt.eidx_filename);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_preproc_param) {
    int argc = 8;
    char* argv[8] = {"build",
                      "aaa.afa", "-o", "aaa.mrf",
                      "--seqwt", "0",
                      "--effnum", "1"};
    MRFBuildCommandLine cmd_line(argc, argv);

    EXPECT_EQ(NO_WEIGHT, cmd_line.opt.msa_analyzer_opt.seq_wt);
    EXPECT_EQ(EXP_ENTROPY, cmd_line.opt.msa_analyzer_opt.eff_num);
}

class MRFStatCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFStatCommandLine_Test, test_parse_param) {
    int argc = 2;
    char* argv[2] = {"stat",
                     "aaa.mrf"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.mrf", cmd_line.opt.mrf_filename);

    EXPECT_EQ(Stat::MODE_PAIR, cmd_line.opt.mode);
    EXPECT_EQ(Stat::CORR_APC, cmd_line.opt.corr);
}

TEST_F(MRFStatCommandLine_Test, test_parse_opt_param) {
    int argc = 6;
    char* argv[6] = {"stat",
                     "aaa.mrf",
                     "--mode", "pos",
                     "--corr", "2"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.mrf", cmd_line.opt.mrf_filename);

    EXPECT_EQ(Stat::MODE_POS, cmd_line.opt.mode);
    EXPECT_EQ(Stat::CORR_NCPS, cmd_line.opt.corr);
}

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
