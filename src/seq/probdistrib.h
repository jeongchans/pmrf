#ifndef _SEQ_PROBDISTRIB_H_
#define _SEQ_PROBDISTRIB_H_

#include "util/numeric.h"

class EmitProbEstimator {
  public:
    virtual Float1dArray estimate(const Float1dArray& c) const = 0;
};

// Maximum likelihood estimation of emission probablities
class MLEmitProbEstimator : public EmitProbEstimator {
  public:
    MLEmitProbEstimator() {};

    virtual Float1dArray estimate(const Float1dArray& c) const;
};

// Substitution matrix mixture estimation of emission probabilities
// Biological sequence analysis, p117-119
class SMMEmitProbEstimator : public EmitProbEstimator {
  public:
    SMMEmitProbEstimator(const Float1dArray& bgfreq, const Float2dArray& scoremat, const double& admix=10.0);

    virtual Float1dArray estimate(const Float1dArray& c) const;

  private:
    Float2dArray cond_prob;    // conditional probability
    const double admix;
};

#endif
