#include "probdistrib.h"

#include "subsmat.h"

Float1dArray MLEmitProbEstimator::estimate(const Float1dArray& c) const {
    Float1dArray p(c.size());
    p = c;
    scale(p);
    return p;
}

SMMEmitProbEstimator::SMMEmitProbEstimator(const Float1dArray& bgfreq, const Float2dArray& scoremat, const double& admix) : admix(admix) {
    int n = bgfreq.size();
    TargetProbEstimatorGivenBG estimator(bgfreq);
    Float2dArray target_prob(n, n);
    target_prob = estimator.probify(scoremat).second;
    cond_prob.resize(n, n);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cond_prob(j, i) = target_prob(i, j) / bgfreq(i);
}

Float1dArray SMMEmitProbEstimator::estimate(const Float1dArray& c) const {
    int n = c.size();
    Float1dArray f(n);     // maximum likelihood probablity
    f = c;
    scale(f);
    Float1dArray psc(n);   // pseudocount
    for (int i = 0; i < n; ++i)
        psc(i) = sum(f * cond_prob.row(i).transpose());
    Float1dArray p(n);     // probability
    double alpha = sum(c) - 1.;
    if (alpha < 0) alpha = 0.;
    p = alpha * f + admix * psc;
    scale(p);
    return p;
}
