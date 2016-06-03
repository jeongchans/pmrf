#ifndef _LBFGS_H_
#define _LBFGS_H_

#include <cstdlib>

#include "lbfgs/lbfgs.h"

class LBFGS {

  public:

    class Parameter {
      public:
        Parameter() : n(0), x(NULL) {};
        Parameter(const int& n) : n(n), x(NULL), fx(0.) { malloc_x(); init_param(); };
        virtual ~Parameter() { if (x != NULL) lbfgs_free(x); }

        int n;
        lbfgsfloatval_t *x;
        lbfgsfloatval_t fx;
        lbfgs_parameter_t opt_param;

      protected:
        void malloc_x() { if (x != NULL) lbfgs_free(x); x = lbfgs_malloc(n); }
        void init_param() { lbfgs_parameter_init(&opt_param); }
    };

    class ObjectiveFunction {
      public:
        virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) = 0;
        int progress(const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls);
    };

    class Optimizer {
      public:
        int optimize(Parameter *param, ObjectiveFunction *obj_func);

      protected:
        static lbfgsfloatval_t _evaluate(void *instance, const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step)
            { return reinterpret_cast<ObjectiveFunction*>(instance)->evaluate(x, g, n, step); }
        static int _progress(void *instance, const lbfgsfloatval_t *x, const lbfgsfloatval_t *g, const lbfgsfloatval_t fx, const lbfgsfloatval_t xnorm, const lbfgsfloatval_t gnorm, const lbfgsfloatval_t step, int n, int k, int ls)
            { return reinterpret_cast<ObjectiveFunction*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls); }
    };

};

#endif
