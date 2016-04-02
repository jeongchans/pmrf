#ifndef _PARAMETERIZE_H_
#define _PARAMETERIZE_H_

#include "lbfgs.h"

#include "seq/trace.h"

#include "mrf.h"
#include "msaanalyze.h"

using std::string;

namespace RegulMethod {
    enum RegulMethod { NONE, L2, PROFILE };
}

class MRFParameterizer {
  public:

    class Parameter : public LBFGS::Parameter {
      public:

        class Option {
          public:
            Option(const int& linesearch=LBFGS_LINESEARCH_BACKTRACKING, const double& delta=1e-4,const int& past=2, const int& max_iterations=500)
            : linesearch(linesearch), delta(delta), past(past), max_iterations(max_iterations) {};

            int linesearch;
            double delta;
            int past;
            int max_iterations;
        };

        Parameter(const MRF& model, const Option& opt);

        int get_nidx(const int& i, const char& p) const;
        inline int get_eidx(const int& i, const int& j, const char& p, const char& q) const;
        inline int get_eidx(const int& i, const int& j, const int& xi, const char& q) const;
        inline int get_eidx(const int& i, const int& j, const char& p, const int& xj) const;
        inline int get_eidx(const int& i, const int& j, const int& xi, const int& xj) const;

        const Alphabet& abc;
        int num_var;
        int length;
        std::unordered_map<EdgeIndex, int> eidx;
        int n_node;
        int n_edge;
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
            Option(const double& lambda1=0.01, const double& lambda2=0.2, const bool& sc=true) 
            : lambda1(lambda1), lambda2(lambda2), sc(sc) {};

            double lambda1;
            double lambda2;
            bool sc;
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
            Option(const RegulMethod::RegulMethod& regul=RegulMethod::RegulMethod::L2, const double& gap_prob=0.14) : regul(regul), gap_prob(gap_prob) {};

            RegulMethod::RegulMethod regul;
            L2Regularization::Option l2_opt;

            double gap_prob;
        };

        ObjectiveFunction(const TraceVector& traces, Parameter& param, Option& opt, const MSAAnalyzer& msa_analyzer);

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

      private:
        const TraceVector& traces;
        const Parameter& param;
        const Option& opt;

        MatrixXi data;
        VectorXf seq_weight;

        const MSAAnalyzer msa_analyzer;

        L2Regularization l2_func;

        MatrixXf calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq);
        VectorXf logsumexp(const MatrixXf& b);
        VectorXf calc_logz(const MatrixXf& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const double& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const double& sw);

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

#endif
