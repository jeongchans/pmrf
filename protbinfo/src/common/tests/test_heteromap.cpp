#include <gtest/gtest.h>

#include "heteromap.h"

class HeteroMapTest : public testing::Test {
  protected:
    string key;
};

TEST_F(HeteroMapTest, test_float) {
    HeteroMap map;

    key = "k1";
    FloatType val = 1.15;
    map.insert(key, val);

    FloatType r_val = *(FloatType*) map[key].ptr.get();
    EXPECT_EQ(val, r_val);
}

TEST_F(HeteroMapTest, test_Float1dArray) {
    HeteroMap map;

    key = "k1";
    Float1dArray val(3);
    val = 1, 2, 3;
    map.insert(key, val);

    Float1dArray r_val = *(Float1dArray*) map[key].ptr.get();
    ASSERT_EQ(val.size(), r_val.size());
    EXPECT_TRUE(all(val == r_val));
    val(1) = 20;
    EXPECT_NE(val(1), r_val(1));
}

TEST_F(HeteroMapTest, test_Float2dArray) {
    HeteroMap map;

    key = "k1";
    Float2dArray val(2, 3);
    val = 1, 2, 3,
          4, 5, 6;
    map.insert(key, val);

    Float2dArray r_val = *(Float2dArray*) map[key].ptr.get();
    ASSERT_EQ(val.rows(), r_val.rows());
    ASSERT_EQ(val.cols(), r_val.cols());
    EXPECT_TRUE(all(val == r_val));
}
