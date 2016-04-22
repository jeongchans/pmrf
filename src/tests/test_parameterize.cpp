#include <gtest/gtest.h>

#include "parameterize.h"

using std::string;

AminoAcid abc("-", false, true);

MSAAnalyzer::Option msa_analyzer_opt;
MSAAnalyzer msa_analyzer(msa_analyzer_opt, abc);

class MRFParameterizer_SymmParameter_Test : public testing::Test {
  public:
    MRFParameterizer_SymmParameter_Test() : length(16), mrf(length, abc) {};

    size_t length;
    MRF mrf;
    MRFParameterizer::Parameter::Option optim_opt;
};

TEST_F(MRFParameterizer_SymmParameter_Test, test_get_nidx) {
    MRFParameterizer::SymmParameter param(mrf, optim_opt);

    int offset = param.v_offset[3];
    EXPECT_EQ(0, offset % 21);
}

TEST_F(MRFParameterizer_SymmParameter_Test, test_get_eidx) {
    MRFParameterizer::SymmParameter param(mrf, optim_opt);
    int i = 2;
    int j = 3;

    int offset = param.w_offset[EdgeIndex(i, j)];
    EXPECT_EQ(0, (offset - param.w_beg_pos()) % (21 * 21));
}

class MRFParameterizer_AsymParameter_Test : public testing::Test {
  public:
    MRFParameterizer_AsymParameter_Test() : length(16), mrf(length, abc) {};

    size_t length;
    MRF mrf;
    MRFParameterizer::Parameter::Option optim_opt;
};

TEST_F(MRFParameterizer_AsymParameter_Test, test_get_eidx) {
    MRFParameterizer::AsymParameter param(mrf, optim_opt);
    int i = 2;
    int j = 3;

    int offset = param.w_offset[EdgeIndex(i, j)];
    EXPECT_EQ(0, (offset - param.w_beg_pos()) % (21 * 21));
}

class MRFParameterizer_RegularizationFunction_Test : public testing::Test {
  protected:
    MRFParameterizer_RegularizationFunction_Test() : length(16), mrf(length, abc), param(mrf, optim_opt), g(NULL) {};

    virtual void SetUp() {
        g = lbfgs_malloc(param.n);
        for (int i = 0; i < param.n; ++i) {
            param.x[i] = exp(randn());
            g[i] = 0.;
        }
        traces.push_back(Trace("MMMMMMMMMMMMMMMMM", "PPDQEFLRGARVQLGDA"));
//        traces.push_back(Trace("UUMMOEEMMIIIMMMMUUUU", "DQHGNRIVHLQ"));
        traces.push_back(Trace("DDMMDDDMMIIIMMMMDDDD", "DQHGNRIVHLQ"));
    }

    virtual void TearDown() {
        lbfgs_free(g);
    }

    lbfgsfloatval_t abs_sum(lbfgsfloatval_t *g) {
        double s = 0.;
        for (int i = 0; i < param.n; ++i) s+= fabs(g[i]);
        return s;
    }

    MRFParameterizer::Parameter::Option optim_opt;
    size_t length;
    MRF mrf;
    MRFParameterizer::SymmParameter param;
    lbfgsfloatval_t *g;
    TraceVector traces;
};

TEST_F(MRFParameterizer_RegularizationFunction_Test, test_l2_regularization) {
    MRFParameterizer::L2Regularization::Option opt(0.01, 0.2);
    MRFParameterizer::L2Regularization regul_func(param, opt);

    lbfgsfloatval_t fx = regul_func.evaluate(param.x, g, param.n, 0);
    EXPECT_TRUE(fx > 0);
    EXPECT_NE(0., abs_sum(g));
}

class MRFParameterizer_Pseudolikelihood_Test : public testing::Test {
  protected:
    MRFParameterizer_Pseudolikelihood_Test() : length(16), mrf(length, abc), param(mrf, optim_opt), g(NULL) {};

    virtual void SetUp() {
        traces.push_back(Trace("MMMMMMMMMMMMMMMMM", "PPDQEFLRGARVQLGDA"));
//        traces.push_back(Trace("UUMMOEEMMIIIMMMMUUUU", "DQHGNRIVHLQ"));
        traces.push_back(Trace("DDMMDDDMMIIIMMMMDDDD", "DQHGNRIVHLQ"));
        seq_weight.resize(2);
        seq_weight << 0.5, 0.5;

        g = lbfgs_malloc(param.n);
        for (int i = 0; i < param.n; ++i) {
            param.x[i] = exp(randn());
            g[i] = 0.;
        }
    }

    virtual void TearDown() {
        lbfgs_free(g);
    }

    MRFParameterizer::Parameter::Option optim_opt;
    size_t length;
    MRF mrf;
    MRFParameterizer::SymmParameter param;
    lbfgsfloatval_t *g;

    TraceVector traces;
    VectorXf seq_weight;
};

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_calc_logpot) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    string seq = traces[0].get_matched_aseq();

    MatrixXl logpot = pll.calc_logpot(param.x, 0, seq);
    EXPECT_EQ(param.num_var, logpot.rows());
    EXPECT_EQ(param.num_node, logpot.cols());
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_calc_logz) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    MatrixXl logpot = randn_matrix(param.num_var, param.num_node).unaryExpr(&exp);

    VectorXl logz = pll.calc_logz(logpot);
    EXPECT_EQ(param.num_node, logz.size());
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_update_obj_score) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);
    MatrixXl logpot = randn_matrix(param.num_var, param.num_node).unaryExpr(&exp);
    VectorXl logz = pll.calc_logz(logpot);
    lbfgsfloatval_t fx = 0.;

    pll.update_obj_score(fx, logpot, logz, 0, seq, sw);
    EXPECT_NE(0., fx);
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_update_gradient) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);
    MatrixXl logpot = randn_matrix(param.num_var, param.num_node).unaryExpr(&exp);
    VectorXl logz = pll.calc_logz(logpot);

    pll.update_gradient(param.x, g, logpot, logz, 0, seq, sw);
    double s = 0.;
    for (int i = 0; i < param.n; ++i) s+= fabs(g[i]);
    EXPECT_NE(0., s);
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_evaluate) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    lbfgsfloatval_t fx = pll.evaluate(param.x, g, param.n, 0);
    EXPECT_NE(0.0, fx);
}

class MRFParameterizer_Test : public testing::Test {
  protected:
    virtual void SetUp() {
        length = 16;
        for (int i = 0; i < 10; ++i) {
            string seq = generate_rand_trace(length);
            traces.push_back(Trace("MMMMMMMMMMMMMMMM", seq));
        }
        for (int i = 0; i < 10; ++i) {
            string seq = generate_rand_trace(length, true);
            ambiguous_traces.push_back(Trace("MMMMMMMMMMMMMMMM", seq));
        }
    }

    string generate_rand_trace(const size_t& length, const bool& ambiguous=false) {
        string letters = "ACDEFGHIKLMNPQRSTVWY";
        if (ambiguous) letters += "BJZOUX";
        size_t n = letters.size();
        string seq = "";
        for (size_t i = 0; i < length; ++i)
            seq += letters[rand() % n];
        return seq;
    }

    MRFParameterizer::Parameter::Option optim_opt;
    size_t length;
    TraceVector traces;
    TraceVector ambiguous_traces;
};

TEST_F(MRFParameterizer_Test, test_update_model) {
    MRF mrf(length, abc);
    MRFParameterizer::SymmParameter param(mrf, optim_opt);
    MRFParameterizer parameterizer(msa_analyzer);
    int i = 2;

    mrf.get_node(i).set_weight(VectorXf::Zero(param.num_var));
    ASSERT_EQ(0., mrf.get_node(i).get_weight().squaredNorm());

    for (int k = 0; k < param.n; ++k) param.x[k] = exp(randn());
    parameterizer.update_model(mrf, param);
    EXPECT_NE(0., mrf.get_node(i).get_weight().squaredNorm());
}

TEST_F(MRFParameterizer_Test, test_parameterize) {
    MRF mrf(length, abc);
    MRFParameterizer parameterizer(msa_analyzer);

    parameterizer.parameterize(mrf, traces);
    EXPECT_NE(0., mrf.get_node(2).get_weight().squaredNorm());
    EXPECT_NE(0., mrf.get_edge(2, 5).get_weight().squaredNorm());
}

TEST_F(MRFParameterizer_Test, test_parameterize_ambiguous) {
    MRF mrf(length, abc);
    MRFParameterizer parameterizer(msa_analyzer);

    parameterizer.parameterize(mrf, ambiguous_traces);
    EXPECT_NE(0., mrf.get_node(2).get_weight().squaredNorm());
    EXPECT_NE(0., mrf.get_edge(2, 5).get_weight().squaredNorm());
}
