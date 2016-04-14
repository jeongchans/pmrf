#ifndef _PARAMETERIZE_H_
#define _PARAMETERIZE_H_

#include "lbfgs.h"

#include "seq/trace.h"

#include "mrf.h"
#include "msaanalyze.h"

using std::string;

namespace RegulMethod {
    enum RegulMethod { REGUL_NONE, REGUL_L2, REGUL_PROFILE };
}

class MRFParameterizer {
  public:

    class Parameter : public LBFGS::Parameter {
      public:

        class Option {
          public:
            Option()
            : linesearch(LBFGS_LINESEARCH_BACKTRACKING),
              delta(1e-4), 
              past(2), 
              max_iterations(500) {};

            int linesearch;
            float delta;
            int past;
            int max_iterations;
        };

        Parameter(const MRF& model, const Option& opt);

        int get_nidx(const int& i, const char& p) const;

        int get_eidx(const int& i, const int& j, const char& p, const char& q) const;
        int get_eidx(const int& i, const int& j, const int& xi, const char& q) const;
        int get_eidx(const int& i, const int& j, const char& p, const int& xj) const;
        int get_eidx(const int& i, const int& j, const int& xi, const int& xj) const;
        int get_eidx_edge(const int& ei, const char& p, const char& q) const;
        int get_eidx_edge(const int& ei, const int& xi, const char& q) const;
        int get_eidx_edge(const int& ei, const char& p, const int& xj) const;
        int get_eidx_edge(const int& ei, const int& xi, const int& xj) const;

        int nidx_beg() const { return 0; }
        int nidx_end() const { return n_node; }
        int eidx_beg() const { return n_node; }
        int eidx_end() const { return n; }

        const Alphabet& abc;
        int length;
        int num_var;
        int n_node;
        int num_edge;
        int num_var2;
        int n_edge;

        std::unordered_map<EdgeIndex, int> eidx;
        EdgeIndexVector edge_idxs;

    };

    // Regularization

    class RegularizationFunction {
      public:
        virtual void regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) = 0;
    };

    // L2 regularization (using zero-mean Gaussian prior)

    class L2Regularization : public RegularizationFunction {
      public:

        class Option {
          public:
            Option(const float& lambda1, const float& lambda2)
            : lambda1(lambda1), lambda2(lambda2) {};

            float lambda1;
            float lambda2;
        };

        L2Regularization(Parameter& param, Option& opt) : param(param), opt(opt) {};

        virtual void regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);

      private:
        const Parameter& param;
        const Option& opt;

        void regularize_node(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);
        void regularize_edge(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);
    };

    // Objective function

    class ObjectiveFunction : public LBFGS::ObjectiveFunction {
      public:

        class Option {
          public:
            Option()
            : regul(RegulMethod::REGUL_L2), 
              l2_opt(regnode_lambda, regedge_lambda),
              regnode_lambda(0.01), 
              regedge_lambda(0.2),
              regedge_sc_deg(true),
              regedge_sc_neff(true) {};

            RegulMethod::RegulMethod regul;
            L2Regularization::Option l2_opt;

            float regnode_lambda;
            float regedge_lambda;
            bool regedge_sc_deg;
            bool regedge_sc_neff;
        };

        ObjectiveFunction(const TraceVector& traces, Parameter& param, Option& opt, const VectorXf& seq_weight, const float& neff);

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

      private:
        const TraceVector& traces;
        const Parameter& param;
        const Option& opt;

        MatrixXi data;
        const VectorXf& seq_weight;
        const float& neff;

        L2Regularization l2_func;

        MatrixXf calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq);
        VectorXf logsumexp(const MatrixXf& b);
        VectorXf calc_logz(const MatrixXf& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw);

        FRIEND_TEST(MRFParameterizer_ObjectiveFunction_Test, test_calc_logpot);
        FRIEND_TEST(MRFParameterizer_ObjectiveFunction_Test, test_logsumexp);
        FRIEND_TEST(MRFParameterizer_ObjectiveFunction_Test, test_calc_logz);
        FRIEND_TEST(MRFParameterizer_ObjectiveFunction_Test, test_update_obj_score);
        FRIEND_TEST(MRFParameterizer_ObjectiveFunction_Test, test_update_gradient);
    };

  public:
    MRFParameterizer(const MSAAnalyzer& msa_analyzer) 
    : msa_analyzer(msa_analyzer) {};

    int parameterize(MRF& model, const TraceVector& traces);

    ObjectiveFunction::Option opt;
    Parameter::Option optim_opt;

  private:
    const MSAAnalyzer& msa_analyzer;

    void update_model(MRF& model, Parameter& param);
    MatrixXf calc_profile(const TraceVector& traces);

    FRIEND_TEST(MRFParameterizer_Test, test_update_model);
};

inline int MRFParameterizer::Parameter::get_nidx(const int& i, const char& p) const 
    { return i * num_var + abc.get_idx(p); }

inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const char& q) const
    { return get_eidx(i, j, abc.get_idx(p), abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const int& xi, const char& q) const
    { return get_eidx(i, j, xi, abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const int& xj) const
    { return get_eidx(i, j, abc.get_idx(p), xj); }
inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const int& xi, const int& xj) const
    { return get_eidx_edge(eidx.at(EdgeIndex(i, j)), xi, xj); }

inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const char& p, const char& q) const
    { return get_eidx_edge(ei, abc.get_idx(p), abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const int& xi, const char& q) const
    { return get_eidx_edge(ei, xi, abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const char& p, const int& xj) const
    { return get_eidx_edge(ei, abc.get_idx(p), xj); }
inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const int& xi, const int& xj) const
    { return eidx_beg() + ei * num_var2 + xi * num_var + xj; }

#endif
