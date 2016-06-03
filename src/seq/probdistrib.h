#ifndef _SEQ_PROBDISTRIB_H_
#define _SEQ_PROBDISTRIB_H_

#include "util/numeric.h"

class EmitProbEstimator {
  public:
    virtual VectorXf estimate(const VectorXf& c) const = 0;
};

// Maximum likelihood estimation of emission probablities
class MLEmitProbEstimator : public EmitProbEstimator {
  public:
    MLEmitProbEstimator() {};

    virtual VectorXf estimate(const VectorXf& c) const;
};

// Substitution matrix mixture estimation of emission probabilities
// Biological sequence analysis, p117-119
class SMMEmitProbEstimator : public EmitProbEstimator {
  public:
    SMMEmitProbEstimator(const VectorXf& bgfreq, const MatrixXf& scoremat, const double& admix=10.0);

    virtual VectorXf estimate(const VectorXf& c) const;

  private:
    MatrixXf cond_prob;    // conditional probability
    const double admix;
};

#endif
