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

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_stat) {
    int argc = 3;
    char* argv[3] = {"pmrf", "stat",
                     "aaa.mrf"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(STAT, cmd_line.subcmd);
}

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_infer) {
    int argc = 4;
    char* argv[4] = {"pmrf", "infer",
                     "aaa.mrf", "aaa.fa"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(INFER, cmd_line.subcmd);
}

TEST_F(MRFMainCommandLine_Test, test_parse_subcmd_show) {
    int argc = 3;
    char* argv[3] = {"pmrf", "show",
                     "aaa.mrf"};
    MRFMainCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(SHOW, cmd_line.subcmd);
}

class MRFBuildCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFBuildCommandLine_Test, test_parse_param_default) {
    int argc = 4;
    char* argv[4] = {"build", "example.afa",
                     "-o", "example.mrf"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("example.afa", cmd_line.opt.msa_filename);
    EXPECT_EQ("example.mrf", cmd_line.opt.out_filename);

    EXPECT_EQ(AFASTA, cmd_line.opt.msa_fmt);

    EXPECT_EQ(MSAProcOption::SW_CLSTR, cmd_line.opt.msa_analyzer_opt.seq_wt);
//    EXPECT_EQ(MSAProcOption::NEFF_CLSTR, cmd_line.opt.msa_analyzer_opt.eff_num);
    EXPECT_FLOAT_EQ(0.8, cmd_line.opt.msa_analyzer_opt.clstr_maxidt);

    EXPECT_EQ(true, cmd_line.opt.parameterizer_opt.asymmetric);
    EXPECT_EQ(RegulMethod::REGUL_L2, cmd_line.opt.parameterizer_opt.regul);
    EXPECT_FLOAT_EQ(UNDETERMINED_F, cmd_line.opt.parameterizer_opt.regnode_lambda);
    EXPECT_FLOAT_EQ(UNDETERMINED_F, cmd_line.opt.parameterizer_opt.regedge_lambda);
//    EXPECT_EQ(true, cmd_line.opt.parameterizer_opt.regedge_sc_deg);
//    EXPECT_EQ(true, cmd_line.opt.parameterizer_opt.regedge_sc_neff);

    EXPECT_EQ(100, cmd_line.opt.optim_opt.corr);
    EXPECT_FLOAT_EQ(1e-5, cmd_line.opt.optim_opt.epsilon);
    EXPECT_FLOAT_EQ(1e-7, cmd_line.opt.optim_opt.delta);
    EXPECT_EQ(500, cmd_line.opt.optim_opt.max_iter);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_input_param) {
    int argc = 8;
    char* argv[8] = {"build", "example.a3m",
                     "-o", "example.mrf",
                     "--msa", "a3m",
                     "--edge", "edges.txt"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("example.a3m", cmd_line.opt.msa_filename);
    EXPECT_EQ("example.mrf", cmd_line.opt.out_filename);
    EXPECT_EQ(A3M, cmd_line.opt.msa_fmt);
    EXPECT_EQ("edges.txt", cmd_line.opt.eidx_filename);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_preproc_param) {
    int argc1 = 4;
    char* argv1[4] = {"build", "example.afa",
                      "--seqwt", "no"};
    MRFBuildCommandLine cmd_line1(argc1, argv1);
    ASSERT_TRUE(cmd_line1.is_valid());
    EXPECT_EQ(MSAProcOption::SW_NO, cmd_line1.opt.msa_analyzer_opt.seq_wt);

    int argc2 = 4;
    char* argv2[4] = {"build", "example.afa",
                      "--clstr-maxidt", "0.4"};
    MRFBuildCommandLine cmd_line2(argc2, argv2);
    ASSERT_TRUE(cmd_line2.is_valid());
    EXPECT_FLOAT_EQ(0.4, cmd_line2.opt.msa_analyzer_opt.clstr_maxidt);
}

//TEST_F(MRFBuildCommandLine_Test, test_parse_neff_clstr_param) {
//    int argc = 8;
//    char* argv[8] = {"build",
//                      "aaa.afa", "-o", "aaa.mrf",
//                      "--neff", "clstr",
//                      "--clstr-maxidt", "0.2"};
//    MRFBuildCommandLine cmd_line(argc, argv);
//
//    EXPECT_EQ(MSAProcOption::NEFF_CLSTR, cmd_line.opt.msa_analyzer_opt.eff_num);
//    EXPECT_FLOAT_EQ(0.2, cmd_line.opt.msa_analyzer_opt.clstr_maxidt);
//}

TEST_F(MRFBuildCommandLine_Test, test_parse_learning_param) {
    int argc = 3;
    char* argv[3] = {"build", "example.afa",
                      "--symmetric"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(false, cmd_line.opt.parameterizer_opt.asymmetric);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_regul_param) {
    int argc1 = 4;
    char* argv1[4] = {"build", "example.afa",
                      "--regul", "no"};
    MRFBuildCommandLine cmd_line1(argc1, argv1);
    ASSERT_TRUE(cmd_line1.is_valid());
    EXPECT_EQ(RegulMethod::REGUL_NONE, cmd_line1.opt.parameterizer_opt.regul);

    int argc2 = 6;
    char* argv2[6] = {"build", "example.afa",
                      "--regul-vl", "0.2",
                      "--regul-wl", "0.02"};
    MRFBuildCommandLine cmd_line2(argc2, argv2);
    ASSERT_TRUE(cmd_line2.is_valid());
    EXPECT_FLOAT_EQ(0.2, cmd_line2.opt.parameterizer_opt.regnode_lambda);
    EXPECT_FLOAT_EQ(0.02, cmd_line2.opt.parameterizer_opt.regedge_lambda);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_lbfgs_param) {
    int argc = 10;
    char* argv[10] = {"build", "example.afa",
                     "--lbfgs-corr", "5",
                     "--lbfgs-epsilon", "1e-3",
                     "--lbfgs-delta", "1e-2",
                     "--lbfgs-maxiter", "100"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ(5, cmd_line.opt.optim_opt.corr);
    EXPECT_FLOAT_EQ(1e-3, cmd_line.opt.optim_opt.epsilon);
    EXPECT_FLOAT_EQ(1e-2, cmd_line.opt.optim_opt.delta);
    EXPECT_EQ(100, cmd_line.opt.optim_opt.max_iter);
}

TEST_F(MRFBuildCommandLine_Test, test_parse_experimental_param) {
    int argc = 8;
    char* argv[8] = {"build", "example.afa",
                     "--regw-lambda-max", "0.4",
                     "--regw-lambda-min", "0.1",
                     "--regw-lambda-sc", "2.0"};
    MRFBuildCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_FLOAT_EQ(0.4, cmd_line.opt.parameterizer_opt.regedge_lambda_max);
    EXPECT_FLOAT_EQ(0.1, cmd_line.opt.parameterizer_opt.regedge_lambda_min);
    EXPECT_FLOAT_EQ(2.0, cmd_line.opt.parameterizer_opt.regedge_lambda_sc);
}

class MRFStatCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFStatCommandLine_Test, test_parse_param_default) {
    int argc = 2;
    char* argv[2] = {"stat", "example.mrf"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("example.mrf", cmd_line.opt.mrf_filename);
    EXPECT_EQ(Stat::MODE_PAIR, cmd_line.opt.mode);
    EXPECT_EQ(Stat::CORR_APC, cmd_line.opt.corr);
}

TEST_F(MRFStatCommandLine_Test, test_parse_pos_param) {
    int argc = 4;
    char* argv[4] = {"stat", "example.mrf",
                     "--mode", "pos"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("example.mrf", cmd_line.opt.mrf_filename);
    EXPECT_EQ(Stat::MODE_POS, cmd_line.opt.mode);
}

TEST_F(MRFStatCommandLine_Test, test_parse_pair_param) {
    int argc = 6;
    char* argv[6] = {"stat", "example.mrf",
                     "--mode", "pair",
                     "--corr", "ncps"};
    MRFStatCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("example.mrf", cmd_line.opt.mrf_filename);
    EXPECT_EQ(Stat::MODE_PAIR, cmd_line.opt.mode);
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

class MRFShowCommandLine_Test : public testing::Test {
  protected:
};

TEST_F(MRFShowCommandLine_Test, test_parse_param) {
    int argc = 2;
    char* argv[2] = {"show",
                     "aaa.mrf"};
    MRFShowCommandLine cmd_line(argc, argv);
    ASSERT_TRUE(cmd_line.is_valid());
    EXPECT_EQ("aaa.mrf", cmd_line.opt.mrf_filename);
}
