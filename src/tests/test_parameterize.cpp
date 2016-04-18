#include <gtest/gtest.h>

#include "parameterize.h"

using std::string;

AminoAcid abc("-", false, true);

MSAAnalyzer::Option msa_analyzer_opt;
MSAAnalyzer msa_analyzer(msa_analyzer_opt, abc);

class MRFParameterizer_Parameter_Test : public testing::Test {
  public:
    MRFParameterizer_Parameter_Test() : length(16), mrf(length, abc) {};

    size_t length;
    MRF mrf;
    MRFParameterizer::Parameter::Option optim_opt;
};

TEST_F(MRFParameterizer_Parameter_Test, test_get_nidx) {
    MRFParameterizer::Parameter param(mrf, optim_opt);

    int nidx = param.get_nidx(3, 'C');
    EXPECT_EQ(3 * 21 + 1, nidx);
}

TEST_F(MRFParameterizer_Parameter_Test, test_get_eidx) {
    MRFParameterizer::Parameter param(mrf, optim_opt);
    int i = 2;
    int j = 3;
    char p = 'E';
    char q = 'C';

    int eidx = param.get_eidx(i, j, p, q);
    int expected = length * 21;
    EdgeIndexVector edge_idxs = mrf.get_edge_idxs();
    for (EdgeIndexVector::iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        if (pos->idx1 == (size_t) i && pos->idx2 == (size_t) j) break;
        else expected += 21 * 21;
    }
    string letters = abc.get_canonical();
    for (size_t k = 0; k < letters.size(); ++k) {
        if (letters[k] == p) break;
        else expected += 21;
    }
    for (size_t k = 0; k < letters.size(); ++k) {
        if (letters[k] == q) break;
        else expected++;
    }

    EXPECT_EQ(expected, eidx);
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
    MRFParameterizer::Parameter param;
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
    MRFParameterizer::Parameter param;
    lbfgsfloatval_t *g;

    TraceVector traces;
    VectorXf seq_weight;
};

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_calc_logpot) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    string seq = traces[0].get_matched_aseq();

    MatrixXf logpot = pll.calc_logpot(param.x, 0, seq);
    EXPECT_EQ(param.num_var, logpot.rows());
    EXPECT_EQ(param.length, logpot.cols());
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_logsumexp) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    MatrixXf x = MatrixXf::Zero(3, 4);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            x(i, j) = 4 * i + j + 1;

    VectorXf y = pll.logsumexp(x);
    ASSERT_TRUE(allclose(y(0), 4.4402));
    ASSERT_TRUE(allclose(y(1), 8.4402));
    ASSERT_TRUE(allclose(y(2), 12.4402));
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_calc_logz) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    MatrixXf logpot = randn_matrix(param.num_var, param.length).unaryExpr(&exp);

    VectorXf logz = pll.calc_logz(logpot);
    EXPECT_EQ(param.length, logz.size());
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_update_obj_score) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);
    MatrixXf logpot = randn_matrix(param.num_var, param.length).unaryExpr(&exp);
    VectorXf logz = pll.calc_logz(logpot);
    lbfgsfloatval_t fx = 0.;

    pll.update_obj_score(fx, logpot, logz, 0, seq, sw);
    EXPECT_NE(0., fx);
}

TEST_F(MRFParameterizer_Pseudolikelihood_Test, test_update_gradient) {
    MRFParameterizer::Pseudolikelihood pll(traces, param, seq_weight, traces.size());

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);
    MatrixXf logpot = randn_matrix(param.num_var, param.length).unaryExpr(&exp);
    VectorXf logz = pll.calc_logz(logpot);

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
    MRFParameterizer::Parameter param(mrf, optim_opt);
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
