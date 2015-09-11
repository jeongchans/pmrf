#include <gtest/gtest.h>

#include "rootfinder.h"

class RootFinderTest : public testing::Test {
  public:
    virtual void SetUp() {
        a = 5.;
        b = 9.;
        c = -2.;
    }

    double a, b, c;
};

TEST_F(RootFinderTest, test_newton_raphson) {
    QuadraticFunction target_func(a, b, c);
    NewtonRaphsonRootFinder solver;
    EXPECT_TRUE(allclose(0.2, solver.find_root(target_func, 1.)));
    EXPECT_TRUE(allclose(-2., solver.find_root(target_func, -3.)));
}
