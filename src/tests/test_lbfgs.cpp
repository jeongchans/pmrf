#include <gtest/gtest.h>

#include "lbfgs.h"

#include "protbinfo/common.h"

class LBFGS_Parameter_Test : public testing::Test {
  protected:
};

TEST_F(LBFGS_Parameter_Test, test_initialize) {
    int n = 10;
    LBFGS::Parameter param(n);
    EXPECT_EQ(n, param.n);
}

class MyObjFunc : public LBFGS::ObjectiveFunction {
  public:
    virtual lbfgsfloatval_t evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
        lbfgsfloatval_t fx = 0.0;
        for (int i = 0; i < n; i += 2) {
            lbfgsfloatval_t t1 = 1.0 - x[i];
            lbfgsfloatval_t t2 = 10.0 * (x[i + 1] - x[i] * x[i]);
            g[i + 1] = 20.0 * t2;
            g[i] = -2.0 * (x[i] * g[i + 1] + t1);
            fx += t1 * t1 + t2 * t2;
        }
        return fx;
    }
};

class LBFGS_Optimizer_Test : public testing::Test {
  protected:
};

TEST_F(LBFGS_Parameter_Test, test_optimize) {
    const int N = 100;
    LBFGS::Optimizer optimizer;
    LBFGS::Parameter param(N);
    MyObjFunc obj_func;
    param.x[0] = 0;
    param.x[1] = 0;

    int ret = optimizer.optimize(&param, &obj_func);
    ASSERT_TRUE(ret >= 0);
    ASSERT_TRUE(param.fx < 1e-10) << param.fx;
    ASSERT_TRUE(allclose(1, param.x[0]));
    ASSERT_TRUE(allclose(1, param.x[1]));
}
