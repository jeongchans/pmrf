#ifndef _ROOTFINDER_H_
#define _ROOTFINDER_H_

#include "common.h"
#include "numeric.h"

class TargetFunction {
  public:
    virtual double fx(const double& x) const = 0;
    virtual pair<double, double> fx_dfx(const double& x) const = 0;
};

// f(x)  = ax^2 + bx + c
// df(x) = 2ax + b
class QuadraticFunction : public TargetFunction {
  public:
    QuadraticFunction(const double& a, const double& b, const double& c)
        : a(a), b(b), c(c) {};

    virtual double fx(const double& x) const {
        return a * x * x + b * x + c;
    }
    virtual pair<double, double> fx_dfx(const double& x) const {
        return make_pair(fx(x), 2 * a * x + b);
    }

  private:
    double a, b, c;
};

class NewtonRaphsonRootFinder {
  public:
    NewtonRaphsonRootFinder() : max_iter(100), atol(1e-15), rtol(1e-15) {};

    double find_root(const TargetFunction& target, const double& guess);

  private:
    int max_iter;
    double atol;
    double rtol;
};

#endif
