#include <gtest/gtest.h>

#include "numeric.h"

TEST(NumericTest, test_raw_data) {
    int rows = 3;
    int cols = 2;
    Float2dArray x = randn(rows, cols);
    std::vector<const FloatType*> r = raw_data(x);
    ASSERT_EQ((size_t) rows, r.size());
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j <cols; ++j)
            EXPECT_EQ(r[i][j], x(i, j));
}

TEST(NumericTest, test_dot) {
    Float2dArray x = randn(3, 2);
    Float2dArray y = randn(2, 4);
    Float2dArray z = dot(x, y);
    ASSERT_EQ(3, z.rows());
    ASSERT_EQ(4, z.cols());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 4; ++j)
            EXPECT_FLOAT_EQ(blitz::sum(x(i, ALL) * y(ALL, j)), z(i, j));

    Float1dArray a = randn(2);
    Float1dArray b = dot(x, a);
    ASSERT_EQ(3, b.size());
    for (int i = 0; i < 3; ++i)
        EXPECT_FLOAT_EQ(blitz::sum(x(i, ALL) * a), b(i));
}

TEST(NumericTest, test_transpose) {
    int rows = 3;
    int cols = 2;
    Float2dArray x = randn(rows, cols);
    Float2dArray y = transpose(x);
    ASSERT_EQ(cols, y.rows());
    ASSERT_EQ(rows, y.cols());
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            EXPECT_EQ(x(j, i), y(i, j));
}

TEST(NumericTest, test_concatenate) {
    int rows1 = 2;
    int rows2 = 1;
    int cols = 2;
    Float2dArray x = randn(rows1, cols);
    Float2dArray y = randn(rows2, cols);
    Float2dArray z = concatenate(x, y);
    ASSERT_EQ(rows1 + rows2, z.rows());
    ASSERT_EQ(cols, z.cols());
    EXPECT_TRUE(all(x == z(Range(blitz::fromStart, rows1 - 1), ALL)));
    EXPECT_TRUE(all(y == z(Range(rows1, blitz::toEnd), ALL)));
}

TEST(NumericTest, test_flatten) {
    int rows = 3;
    int cols = 2;
    Float2dArray x = randn(rows, cols);
    Float1dArray y = flatten(x);
    ASSERT_EQ(rows * cols, y.size());
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            EXPECT_EQ(y(i * cols + j), x(i, j));
}

TEST(NumericTest, test_export_and_import_float1darray) {
    std::string s;
    int size = 3;
    Float1dArray x(size);
    x = 0.1, 0.2, 0.3;

    std::ostringstream oss;
    export_array(x, oss);
    s = "3\n0.1\t0.2\t0.3\n";
    EXPECT_EQ(s, oss.str());

    std::istringstream iss(s);
    Float1dArray y;
    import_array(y, iss);
    ASSERT_EQ(size, y.size());
    EXPECT_TRUE(all(x == y));
}

TEST(NumericTest, test_export_and_import_float2darray) {
    std::string s;
    int rows = 3;
    int cols = 2;
    Float2dArray x(rows, cols);
    x = 0.1, 0.2,
        0.3, 0.4,
        0.5, 0.6;

    std::ostringstream oss;
    export_array(x, oss);
    s = "3 x 2\n"
        "0.1\t0.2\n"
        "0.3\t0.4\n"
        "0.5\t0.6\n";
    EXPECT_EQ(s, oss.str());

    std::istringstream iss(s);
    Float2dArray y;
    import_array(y, iss);
    ASSERT_EQ(rows, y.rows());
    ASSERT_EQ(cols, y.cols());
    EXPECT_TRUE(all(x == y));
}

class RowColReductionTest : public testing::Test {
  protected:
    virtual void SetUp() {
        rows = 2;
        cols = 3;
        x.resize(rows, cols);
        x = 1, 2, 3,
            4, 5, 6;
    }

    Float2dArray x;
    int rows, cols;
};

TEST_F(RowColReductionTest, test_row_sum) {
    Float1dArray s = row_sum(x);
    ASSERT_EQ(rows, s.size());
    EXPECT_EQ(6, s(0));
    EXPECT_EQ(15, s(1));
}

TEST_F(RowColReductionTest, test_col_sum) {
    Float1dArray s = col_sum(x);
    ASSERT_EQ(cols, s.size());
    EXPECT_EQ(5, s(0));
    EXPECT_EQ(7, s(1));
    EXPECT_EQ(9, s(2));
}

TEST_F(RowColReductionTest, test_row_mean) {
    Float1dArray mn = row_mean(x);
    ASSERT_EQ(rows, mn.size());
    EXPECT_EQ(2, mn(0));
    EXPECT_EQ(5, mn(1));
}

TEST_F(RowColReductionTest, test_col_mean) {
    Float1dArray mn = col_mean(x);
    ASSERT_EQ(cols, mn.size());
    EXPECT_EQ(2.5, mn(0));
    EXPECT_EQ(3.5, mn(1));
    EXPECT_EQ(4.5, mn(2));
}
