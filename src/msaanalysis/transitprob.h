#ifndef _PROTBINFO_TRANSITPROB_H_
#define _PROTBINFO_TRANSITPROB_H_

#include "common/numeric.h"
#include "align/trace.h"

class TransitProbEstimator {
  public:
    virtual Float1dArray estimate(const Float1dArray& c, const StateType& t) const = 0;
};

// Maximum likelihood estimation of transition probabilities
class MLTransitProbEstimator : public TransitProbEstimator {
  public:
    MLTransitProbEstimator() {};

    virtual Float1dArray estimate(const Float1dArray& c, const StateType&) const;
};

// Prior based transition probability estimator
class PriorTransitProbEstimator : public TransitProbEstimator {
  public:
    PriorTransitProbEstimator() : admix(1.0) {};
    PriorTransitProbEstimator(const double& admix) : admix(admix) {};
    PriorTransitProbEstimator(const double& mm, const double& mi, const double& md, const double& im, const double& ii, const double& dm, const double& dd, const double& admix);

    virtual Float1dArray estimate(const Float1dArray& c, const StateType& t) const;

  protected:
    Float1dArray prior[NUM_STATE_TYPE];
    double admix;
};

// Empirical maximum discrimination prior based transition probability estimator
// Wistrand and Sonnhammer. J Comput Biol (2004) vol. 11 (1) pp. 181-93
class EMDTransitProbEstimator : public PriorTransitProbEstimator {
  public:
    EMDTransitProbEstimator(const double& admix=1.0);
};

// Transition priors: originally from Graeme Mitchison
// The data is from HMMER-3.0 source code.
class GMTransitProbEstimator : public PriorTransitProbEstimator {
  public:
    GMTransitProbEstimator(const double& admix=1.0);
};

#endif
