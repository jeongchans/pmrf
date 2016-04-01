#include "probdistrib.h"

#include "subsmat.h"

VectorXf MLEmitProbEstimator::estimate(const VectorXf& c) const {
    return c / c.sum();
}

SMMEmitProbEstimator::SMMEmitProbEstimator(const VectorXf& bgfreq, const MatrixXf& scoremat, const double& admix) : admix(admix) {
    int n = bgfreq.size();
    TargetProbEstimatorGivenBG estimator(bgfreq);
    MatrixXf target_prob(n, n);
    target_prob = estimator.probify(scoremat).second;
    cond_prob.resize(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cond_prob(j, i) = target_prob(i, j) / bgfreq(i);
}

VectorXf SMMEmitProbEstimator::estimate(const VectorXf& c) const {
    int n = c.size();
    VectorXf f = c / c.sum();     // maximum likelihood probablity
    VectorXf psc(n);   // pseudocount
    for (int i = 0; i < n; ++i)
        psc(i) = cond_prob.row(i) * f;
    VectorXf p(n);     // probability
    double alpha = c.sum() - 1.;
    if (alpha < 0) alpha = 0.;
    p = alpha * f + admix * psc;
    return p / p.sum();
}
