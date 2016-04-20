#ifndef _PARAMETERIZE_H_
#define _PARAMETERIZE_H_

#include <utility>

#include "lbfgs.h"

#include "seq/trace.h"

#include "mrf.h"
#include "msaanalyze.h"

using std::string;

namespace RegulMethod {
    enum RegulMethod { REGUL_NONE, REGUL_L2, REGUL_PROFILE };
}

typedef LBFGS::ObjectiveFunction* PtrObjFunc;

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
            : corr(10),
              epsilon(1e-5),
              past(2), 
              delta(1e-4), 
              linesearch(LBFGS_LINESEARCH_BACKTRACKING_ARMIJO),
              max_iterations(500) {};

            int corr;
            float epsilon;
            int past;
            float delta;
            int linesearch;
            int max_iterations;
        };

        Parameter(const MRF& model, const Option& opt) : abc(model.get_alphabet()), length(model.get_length()) {};

        int get_nidx(const int& i, const char& p) const;
        int get_nidx(const int& i, const int& xi) const;

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

      protected:
        void init(const MRF& model, const Option& opt);
        void set_opt(const Option& opt);
    };

    class SymmParameter : public Parameter {
      public:
        SymmParameter(const MRF& model, const Option& opt);

        EdgeIndexVector edge_idxs;
    };

    class LocalParameter : public Parameter {
      public:
        LocalParameter(const MRF& model, const Option& opt, const Parameter& p, const int& obs_node);

        int get_nidx(const int& i, const int& xi) const { return xi; }

      private:
        const int& obs_node;
    };

    class AsymParameter : public Parameter {
      public:
        AsymParameter(const MRF& model, const Option& opt);

        void set_x(const LocalParameter& p, const int& obs_node);
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

        MatrixXf calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq);
        VectorXf logsumexp(const MatrixXf& b);
        VectorXf calc_logz(const MatrixXf& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw);

        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_calc_logpot);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_logsumexp);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_calc_logz);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_update_obj_score);
        FRIEND_TEST(MRFParameterizer_Pseudolikelihood_Test, test_update_gradient);
    };

    // Asymmetric pseudolikelihood function

    class AsymPseudolikelihood : public LBFGS::ObjectiveFunction {
      public:
        AsymPseudolikelihood(const TraceVector& traces, const Parameter& param, const VectorXf& seq_weight, const float& neff);

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

        void set_local_param(const Parameter* p) { lp = p; }
        void set_obs_node(const size_t& i) { obs_node = i; }

      private:
        const TraceVector& traces;
        const Parameter& param;
        const Parameter* lp;

        MatrixXi data;
        const VectorXf& seq_weight;
        const float& neff;

        int obs_node;

        VectorXf calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq, const Parameter& param);
        float calc_logz(const VectorXf& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const VectorXf& logpot, const float& logz, const size_t m, const string& seq, const float& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const VectorXf& logpot, const float& logz, const size_t m, const string& seq, const float& sw, const Parameter& param);
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
          regedge_lambda(0.2),
          regedge_sc_deg(true),
          regedge_sc_neff(true),
          asymmetric(true) {};

        RegulMethod::RegulMethod regul;
        float regnode_lambda;
        float regedge_lambda;
        bool regedge_sc_deg;
        bool regedge_sc_neff;
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
    MatrixXf calc_profile(const TraceVector& traces);

    FRIEND_TEST(MRFParameterizer_Test, test_update_model);
};

inline int MRFParameterizer::Parameter::get_nidx(const int& i, const char& p) const 
    { return get_nidx(i, abc.get_idx(p)); }
inline int MRFParameterizer::Parameter::get_nidx(const int& i, const int& xi) const 
    { return i * num_var + xi; }

inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const char& q) const
    { return get_eidx(i, j, abc.get_idx(p), abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const int& xi, const char& q) const
    { return get_eidx(i, j, xi, abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const int& xj) const
    { return get_eidx(i, j, abc.get_idx(p), xj); }
inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const int& xi, const int& xj) const
    { return get_eidx_edge(eidx.at(EdgeIndex((size_t) i, (size_t) j)), xi, xj); }

inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const char& p, const char& q) const
    { return get_eidx_edge(ei, abc.get_idx(p), abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const int& xi, const char& q) const
    { return get_eidx_edge(ei, xi, abc.get_idx(q)); }
inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const char& p, const int& xj) const
    { return get_eidx_edge(ei, abc.get_idx(p), xj); }
inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const int& xi, const int& xj) const
    { return eidx_beg() + ei * num_var2 + xi * num_var + xj; }

#endif
