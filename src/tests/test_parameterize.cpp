#include <gtest/gtest.h>

#include "parameterize.h"

using std::string;

// consider gap as an additional symbol
//Alphabet abc("ACDEFGHIKLMNPQRSTVWY-",
//             "",
//             "BJZOU",
//             "X",
//             "*",
//             "~");

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
        if (pos->idx1 == i && pos->idx2 == j) break;
        else expected += 21 * 21;
    }
    string letters = abc.get_canonical();
    for (int k = 0; k < letters.size(); ++k) {
        if (letters[k] == p) break;
        else expected += 21;
    }
    for (int k = 0; k < letters.size(); ++k) {
        if (letters[k] == q) break;
        else expected++;
    }

    EXPECT_EQ(expected, eidx);
}

class MRFParameterizer_RegularizationFunction_Test : public testing::Test {
  protected:
    MRFParameterizer_RegularizationFunction_Test() : length(16), g(NULL), mrf(length, abc), param(mrf, optim_opt) {};

    virtual void SetUp() {
        g = lbfgs_malloc(param.n);
        for (int i = 0; i < param.n; ++i) {
            param.x[i] = exp(randn());
            g[i] = 0.;
        }
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
};

TEST_F(MRFParameterizer_RegularizationFunction_Test, test_node_l2_regularization) {
    MRFParameterizer::NodeL2Regularization::Option opt;
    MRFParameterizer::NodeL2Regularization regul_func(param, opt);
    lbfgsfloatval_t fx = 0.0;

    regul_func.regularize(param.x, g, fx);
    EXPECT_TRUE(fx > 0);
    EXPECT_NE(0., abs_sum(g));
}

TEST_F(MRFParameterizer_RegularizationFunction_Test, test_edge_l2_regularization) {
    MRFParameterizer::EdgeL2Regularization::Option opt;
    MRFParameterizer::EdgeL2Regularization regul_func(param, opt);
    lbfgsfloatval_t fx = 0.0;

    regul_func.regularize(param.x, g, fx);
    EXPECT_TRUE(fx > 0);
    EXPECT_NE(0., abs_sum(g));
}

class MRFParameterizer_ObjectiveFunction_Test : public testing::Test {
  protected:
    MRFParameterizer_ObjectiveFunction_Test() : length(16), g(NULL), mrf(length, abc), param(mrf, optim_opt) {};

    virtual void SetUp() {
        size_t n = abc.get_canonical_size();
        traces.push_back(Trace("MMMMMMMMMMMMMMMMM", "PPDQEFLRGARVQLGDA"));
        traces.push_back(Trace("UUMMOEEMIIIIMMMMUUUU", "DQHGNRIVHLQ"));
        seq_weight.resize(2);
        seq_weight = 0.5, 0.5;

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
    MRFParameterizer::ObjectiveFunction::Option opt;
    lbfgsfloatval_t *g;

    TraceVector traces;
    Float1dArray seq_weight;
};

TEST_F(MRFParameterizer_ObjectiveFunction_Test, test_calc_logpot) {
    MRFParameterizer::ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);

    Float2dArray logpot = obj_func.calc_logpot(param.x, seq, sw);
    EXPECT_EQ(param.num_var, logpot.rows());
    EXPECT_EQ(param.length, logpot.cols());
}

TEST_F(MRFParameterizer_ObjectiveFunction_Test, test_logsumexp) {
    MRFParameterizer::ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);

    Float2dArray x = zeros(3, 4);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            x(i, j) = 4 * i + j + 1;

    Float1dArray y = obj_func.logsumexp(x);
    ASSERT_TRUE(allclose(y(0), 4.4402));
    ASSERT_TRUE(allclose(y(1), 8.4402));
    ASSERT_TRUE(allclose(y(2), 12.4402));
}

TEST_F(MRFParameterizer_ObjectiveFunction_Test, test_calc_logz) {
    MRFParameterizer::ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);

    Float2dArray logpot = exp(randn(param.num_var, param.length));

    Float1dArray logz = obj_func.calc_logz(logpot);
    EXPECT_EQ(param.length, logz.size());
}

TEST_F(MRFParameterizer_ObjectiveFunction_Test, test_update_obj_score) {
    MRFParameterizer::ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);
    Float2dArray logpot = exp(randn(param.num_var, param.length));
    Float1dArray logz = obj_func.calc_logz(logpot);
    lbfgsfloatval_t fx = 0.;

    obj_func.update_obj_score(fx, logpot, logz, seq, sw);
    EXPECT_NE(0., fx);
}

TEST_F(MRFParameterizer_ObjectiveFunction_Test, test_update_gradient) {
    MRFParameterizer::ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);

    string seq = traces[0].get_matched_aseq();
    FloatType sw = seq_weight(0);
    Float2dArray logpot = exp(randn(param.num_var, param.length));
    Float1dArray logz = obj_func.calc_logz(logpot);

    obj_func.update_gradient(param.x, g, logpot, logz, seq, sw);
    double s = 0.;
    for (int i = 0; i < param.n; ++i) s+= fabs(g[i]);
    EXPECT_NE(0., s);
}

TEST_F(MRFParameterizer_ObjectiveFunction_Test, test_evaluate) {
    MRFParameterizer::ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);

    lbfgsfloatval_t fx = obj_func.evaluate(param.x, g, param.n, 0);
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
        for (int i = 0; i < length; ++i)
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

    mrf.get_node(i).set_weight(zeros(param.num_var));
    ASSERT_EQ(0., sum(pow2(mrf.get_node(i).get_weight())));

    for (int k = 0; k < param.n; ++k) param.x[k] = exp(randn());
    parameterizer.update_model(mrf, param);
    EXPECT_NE(0., sum(pow2(mrf.get_node(i).get_weight())));
}

TEST_F(MRFParameterizer_Test, test_parameterize) {
    MRF mrf(length, abc);
    MRFParameterizer parameterizer(msa_analyzer);

    int ret = parameterizer.parameterize(mrf, traces);
    //ASSERT_TRUE(ret >= 0);
    EXPECT_NE(0., sum(pow2(mrf.get_node(2).get_weight())));
    EXPECT_NE(0., sum(pow2(mrf.get_edge(2, 5).get_weight())));
}

TEST_F(MRFParameterizer_Test, test_parameterize_ambiguous) {
    MRF mrf(length, abc);
    MRFParameterizer parameterizer(msa_analyzer);

    int ret = parameterizer.parameterize(mrf, ambiguous_traces);
    //ASSERT_TRUE(ret >= 0);
    EXPECT_NE(0., sum(pow2(mrf.get_node(2).get_weight())));
    EXPECT_NE(0., sum(pow2(mrf.get_edge(2, 5).get_weight())));
}
