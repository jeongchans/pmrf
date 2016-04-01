#include "rootfinder.h"

#include "util/numeric.h"

double NewtonRaphsonRootFinder::find_root(const TargetFunction& target, const double& guess) {
    double x = guess;
    pair<double, double> p = target.fx_dfx(x);
    double fx = p.first;
    double dfx = p.second;
    double old_x, old_fx;
    int i = 0;
    while (true) {
        ++i;
        if (i > max_iter) 
            throw ConvergeException();
        old_x = x;
        old_fx = fx;
        x = x - fx / dfx;
        p = target.fx_dfx(x);
        fx = p.first;
        dfx = p.second;
        if (fx == 0) break;     // exact solution
        if (fabs(x - old_x) < atol + rtol * x) break;
    }
    return x;
}
