#include <gtest/gtest.h>

#include "mrfio.h"

class EdgeIndexImporter_Test : public testing::Test {
  protected:
    virtual void SetUp() {
        buf =
            "1\t2\n"
            "1\t3\n"
            "8\t10\n"
            "9\t10";
    }

    string buf;
};

TEST_F(EdgeIndexImporter_Test, test_import) {
    std::istringstream is(buf);
    EdgeIndexImporter importer;
    EdgeIndexVector ret = importer.import(is);
    EXPECT_EQ(4, ret.size());
    EXPECT_TRUE(EdgeIndex(0, 1) == ret[0]);
    EXPECT_TRUE(EdgeIndex(0, 2) == ret[1]);
    EXPECT_TRUE(EdgeIndex(7, 9) == ret[2]);
    EXPECT_TRUE(EdgeIndex(8, 9) == ret[3]);
}
