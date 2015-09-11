#include "transitprob.h"

Float1dArray MLTransitProbEstimator::estimate(const Float1dArray& c, const StateType&) const {
    Float1dArray p(c.shape());
    p = c;
    scale(p);
    return p;
}

PriorTransitProbEstimator::PriorTransitProbEstimator(const double& mm, const double& mi, const double& md, const double& im, const double& ii, const double& dm, const double& dd, const double& admix) : admix(admix) {
    for (size_t i = 0; i < NUM_STATE_TYPE; ++i)
        prior[i].resize(3);
    prior[MATCH]  = mm, md, mi;
    prior[DELETE] = dm, dd,  0;
    prior[INSERT] = im,  0, ii;
}

Float1dArray PriorTransitProbEstimator::estimate(const Float1dArray& c, const StateType& t) const {
    Float1dArray p(c.shape());
    double alpha = blitz::sum(c);
    if (alpha < 0) alpha = 0.;
    Float1dArray f = c.copy();
    scale(f);
    p = alpha * f + admix * prior[t];
    scale(p);
    return p;
}

EMDTransitProbEstimator::EMDTransitProbEstimator(const double& admix) : PriorTransitProbEstimator(admix) {
    for (size_t i = 0; i < NUM_STATE_TYPE; ++i)
        prior[i].resize(3);
    prior[MATCH]  = 0.794, 0.005, 0.095;
    prior[DELETE] = 0.278, 0.222, 0.;
    prior[INSERT] = 0.333, 0.,    0.667;
}

GMTransitProbEstimator::GMTransitProbEstimator(const double& admix) : PriorTransitProbEstimator(admix) {
    for (size_t i = 0; i < NUM_STATE_TYPE; ++i)
        prior[i].resize(3);
    prior[MATCH]  = 0.7939, 0.0135, 0.0278;
    prior[DELETE] = 0.9002, 0.5630, 0.;
    prior[INSERT] = 0.1551, 0.,     0.1331;
}
