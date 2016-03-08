#ifndef _PARAMETERIZE_H_
#define _PARAMETERIZE_H_

#include "lbfgs.h"

#include "seq/trace.h"

#include "mrf.h"
#include "msaanalyze.h"

using std::string;

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
        int get_eidx(const int& i, const int& j, const char& p, const char& q) const;

        const Alphabet& abc;
        int num_var;
        int length;
        map<EdgeIndex, int> eidx;
    };

    // Regularization

    class RegularizationFunction {
      public:
        virtual void regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) = 0;
    };

    // Node L2 regularization

    class NodeL2Regularization : public RegularizationFunction {
      public:

        class Option {
          public:
            Option(const double& lambda=0.01) : lambda(lambda) {};

            double lambda;
        };

        NodeL2Regularization(Parameter& param, Option& opt) : param(param), opt(opt) {};

        virtual void regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);

      private:
        const Parameter& param;
        const Option& opt;
    };

    // Edge L2 regularization

    class EdgeL2Regularization : public RegularizationFunction {
      public:

        class Option {
          public:
            Option(const double& lambda=0.2, const bool& sc=true) : lambda(lambda), sc(sc) {};

            double lambda;
            bool sc;
        };

        EdgeL2Regularization(Parameter& param, Option& opt);

        virtual void regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx);

      private:
        const Parameter& param;
        const Option& opt;
        double lambda; 
    };

    // Objective function

    class ObjectiveFunction : public LBFGS::ObjectiveFunction {
      public:

        class Option {
          public:
            Option(const bool& node_l2_regul=true, const bool& edge_l2_regul=true)
            : node_l2_regul(node_l2_regul), edge_l2_regul(edge_l2_regul) {};

            bool node_l2_regul;
            NodeL2Regularization::Option node_l2_opt;
            bool edge_l2_regul;
            EdgeL2Regularization::Option edge_l2_opt;
        };

        ObjectiveFunction(const TraceVector& traces, Parameter& param, Option& opt, const MSAAnalyzer& msa_analyzer);

        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step);

      private:
        const TraceVector& traces;
        const Parameter& param;
        const Option& opt;

        Float1dArray seq_weight;

        const MSAAnalyzer msa_analyzer;

        NodeL2Regularization node_l2_func;
        EdgeL2Regularization edge_l2_func;

        Float2dArray calc_logpot(const lbfgsfloatval_t *x, const string& seq, const double& sw);
        Float1dArray logsumexp(const Float2dArray& b);
        Float1dArray calc_logz(const Float2dArray& logpot);
        void update_obj_score(lbfgsfloatval_t& fx, const Float2dArray& logpot, const Float1dArray& logz, const string& seq, const double& sw);
        void update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const Float2dArray& logpot, const Float1dArray& logz, const string& seq, const double& sw);

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
    Float2dArray calc_profile(const TraceVector& traces);

    FRIEND_TEST(MRFParameterizer_Test, test_update_model);
};

#endif
