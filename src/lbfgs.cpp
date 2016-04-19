#include "lbfgs.h"

#include <iostream>

// ObjectiveFunction

int LBFGS::ObjectiveFunction::progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls) {
//    std::cout << "Iteration " << k;
//    std::cout << ", obj_val = " << fx << std::endl;
    return 0;
}

// Optimizer

int LBFGS::Optimizer::optimize(Parameter *param, ObjectiveFunction *obj_func) {
    int ret = lbfgs(param->n, param->x, &(param->fx), _evaluate, _progress, obj_func, &(param->opt_param));
    if (ret < 0) std::cerr << "L-BFGS error = " << ret << std::endl;
    return ret;
}
