#ifndef _PARAMETERIZE_H_
#define _PARAMETERIZE_H_

#include <utility>

#include "lbfgs.h"

#include "seq/trace.h"

#include "mrf.h"
#include "msaanalyze.h"

#ifdef _TEST_
#include <gtest/gtest_prod.h>
#endif

using std::string;
using std::unordered_map;

namespace RegulMethod {
    enum RegulMethod { REGUL_NONE, REGUL_L2, REGUL_PROFILE };
}

typedef LBFGS::ObjectiveFunction* PtrObjFunc;

typedef Eigen::Matrix<lbfgsfloatval_t, Eigen::Dynamic, 1> VectorXl;
typedef Eigen::Matrix<lbfgsfloatval_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> MatrixXl;
typedef Eigen::Map<VectorXl> MapVectorXl;
typedef Eigen::Map<const VectorXl> MapKVectorXl;
typedef Eigen::Map<MatrixXl> MapMatrixXl;
typedef Eigen::Map<const MatrixXl> MapKMatrixXl;

namespace std {
    template<> struct hash<EdgeIndex> {
        size_t operator()(const EdgeIndex& eidx) const { return eidx.idx1 + eidx.idx2 * 10000; }
    };
}

class MRFParameterizer {
  public:

    class Parameter : public LBFGS::Parameter {
      public:

        class Option {
          public:
            Option()
            : corr(50),
              epsilon(1e-5),
              past(1), 
              delta(1e-5), 
              max_iter(500),
              linesearch(LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE),
              max_linesearch(30) {};

            int corr;
            float epsilon;
            int past;
            float delta;
            int max_iter;
            int linesearch;
            int max_linesearch;;
        };

        Parameter(const Alphabet& abc, const size_t& num_var) 
        : LBFGS::Parameter(), abc(abc), num_var(num_var), num_var2(num_var * num_var) {};

        size_t v_beg_pos() const { return 0; }
        size_t v_end_pos() const { return n_v; }
        size_t w_beg_pos() const { return n_v; }
        size_t w_end_pos() const { return n; }

        const Alphabet& abc;
        const size_t num_var;
        const size_t num_var2;
        size_t num_node;
        size_t num_edge;
        size_t n_v;
        size_t n_w;

        unordered_map<NodeIndex, size_t> v_offset;
        unordered_map<EdgeIndex, size_t> w_offset;

      protected:
        void init(const vector<NodeIndex>& vidxs, const vector<EdgeIndex>& widxs, const Option& opt);
        void set_opt(const Option& opt);
    };

    class SymmParameter : public Parameter {
      public:
        SymmParameter(const MRF& model, const Option& opt);

        EdgeIndexVector edge_idxs;
    };

    class LocalParameter : public Parameter {
      public:
        LocalParameter(const MRF& model, const Option& opt, const Parameter& p, const size_t& obs_node);

        size_t get_obs_node() const { return v_offset.cbegin()->first; }
    };

    class AsymParameter : public Parameter {
      public:
        AsymParameter(const MRF& model, const Option& opt);

        void set_x(const LocalParameter& p);
    };

    // L2 regularization (using zero-mean Gaussian prior)

    class L2Regularization : public LBFGS::ObjectiveFunction {
      public:

        struct Option {
            float lambda1;
            float lambda2;

            Option(const float& lambda1, const float& lambda2) : lambda1(lambda1), lambda2(lambda2) {}
        };

        L2Regularization(Parameter& param, Option& opt) : param(param), opt(opt) {};

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

      private:
        const Parameter& param;
        const Option& opt;

        void regularize_node(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);
        void regularize_edge(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);
    };

    // Pseudolikelihood function

    class Pseudolikelihood : public LBFGS::ObjectiveFunction {
      public:
        Pseudolikelihood(const TraceVector& traces, const SymmParameter& param, const VectorXf& seq_weight, const float& neff);

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

      private:
        const TraceVector& traces;
        const SymmParameter& param;

        MatrixXi data;
        const VectorXf& seq_weight;
        const float& neff;

        MatrixXl calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq);
        VectorXl calc_logz(const MatrixXl& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const MatrixXl& logpot, const VectorXl& logz, const size_t m, const string& seq, const float& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXl& logpot, const VectorXl& logz, const size_t m, const string& seq, const float& sw);

#ifdef _TEST_
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_calc_logpot);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_calc_logz);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_update_obj_score);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_update_gradient);
#endif
    };

    // Asymmetric pseudolikelihood function

    class AsymPseudolikelihood : public LBFGS::ObjectiveFunction {
      public:
        AsymPseudolikelihood(const TraceVector& traces, const Parameter& param, const VectorXf& seq_weight, const float& neff);

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

        void set_local_param(const LocalParameter* p) { lp = p; }

      private:
        const TraceVector& traces;
        const Parameter& param;
        const LocalParameter* lp;

        MatrixXi data;
        vector<string> aseq;
        const VectorXf& seq_weight;
        const float& neff;

        VectorXl calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq);
        lbfgsfloatval_t calc_logz(const VectorXl& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const VectorXl& logpot, const lbfgsfloatval_t& logz, const size_t m, const string& seq, const float& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const VectorXl& logpot, const lbfgsfloatval_t& logz, const size_t m, const string& seq, const float& sw);
    };

    // Objective function

    class ObjectiveFunction : public LBFGS::ObjectiveFunction {
      public:
        ObjectiveFunction(const vector<PtrObjFunc>& funcs, const Parameter& param) 
        : funcs(funcs), param(param) {};

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

      private:
        const vector<PtrObjFunc>& funcs;
        const Parameter& param;
    };

  public:

    class Option {
      public:
        Option()
        : regul(RegulMethod::REGUL_L2), 
          l2_opt(regnode_lambda, regedge_lambda),
          regnode_lambda(0.01), 
          regedge_lambda(0.01),
          asymmetric(true) {};

        RegulMethod::RegulMethod regul;
        float regnode_lambda;
        float regedge_lambda;

        bool asymmetric;

        L2Regularization::Option l2_opt;
    };

    MRFParameterizer(const MSAAnalyzer& msa_analyzer) 
    : msa_analyzer(msa_analyzer) {};

    int parameterize(MRF& model, const TraceVector& traces);

    Option opt;
    Parameter::Option optim_opt;

  private:
    const MSAAnalyzer& msa_analyzer;

    void update_model(MRF& model, const SymmParameter& param);
    void update_model(MRF& model, const AsymParameter& param);
    void adjust_gauge(SymmParameter& param);
    void adjust_gauge(AsymParameter& param);
    void get_reg_lambda(float& regnode_lambda, float& regedge_lambda, const float& avg_deg, const float& neff);

#ifdef _TEST_
    FRIEND_TEST(MRFParameterizer_Test, test_update_model);
#endif
};

#endif
