#include "lbfgs.h"

#include <iostream>

// ObjectiveFunction

int LBFGS::ObjectiveFunction::progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
#ifdef _DEBUG_
    lbfgsfloatval_t opt_cond;
    if (xnorm < 1.) opt_cond = gnorm;
    else opt_cond = gnorm / xnorm;
    std::clog << "[L-BFGS progress]"
              << "  iteration = " << k
              << ", obj_val = " << fx
              << ", xnorm = " << xnorm
              << ", gnorm = " << gnorm
              << ", opt_cond = " << opt_cond
              << ", step = " << step
              << ", ls = " << ls 
              << std::endl;
#endif
    return 0;
}

// Optimizer

int LBFGS::Optimizer::optimize(Parameter *param, ObjectiveFunction *obj_func) {
    int ret = lbfgs(param->n, param->x, &(param->fx), _evaluate, _progress, obj_func, &(param->opt_param));
#ifdef _DEBUG_
    if (ret < 0) {
        std::clog << "[L-BFGS error]"
                  << "  ret = " << ret
                  << std::endl;
    }
#endif
    return ret;
}
