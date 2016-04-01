#include <gtest/gtest.h>

#include "numeric.h"

TEST(NumericTest, test_randn_vector) {
    VectorXf x = randn_vector(3);
    ASSERT_EQ(x.size(), 3);
}

TEST(NumericTest, test_randn_matrix) {
    MatrixXf x = randn_matrix(3, 4);
    ASSERT_EQ(x.rows(), 3);
    ASSERT_EQ(x.cols(), 4);
}

//TEST(NumericTest, test_randn) {
//    Float1dArray x = randn(3);
//    ASSERT_EQ(x.size(), 3);
//}

//TEST(NumericTest, test_transpose) {
//    int rows = 3;
//    int cols = 2;
//    Float2dArray x = randn(rows, cols);
//    Float2dArray y = transpose(x);
//    ASSERT_EQ(cols, y.rows());
//    ASSERT_EQ(rows, y.cols());
//    for (int i = 0; i < rows; ++i)
//        for (int j = 0; j < cols; ++j)
//            EXPECT_EQ(x(i, j), y(j, i));
//}
