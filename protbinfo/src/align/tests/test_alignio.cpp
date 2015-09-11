#include <gtest/gtest.h>

#include "alignio.h"

class FastaAlignmentExporterTest : public testing::Test {
  protected:
    virtual void SetUp() {
        s = "";
        oss.flush();
    }

    FastaAlignmentExporter exporter;
    std::string s;
    std::ostringstream oss;
};

TEST_F(FastaAlignmentExporterTest, test_export_pos_idx) {
    std::vector<int> pos_idx;
    pos_idx.push_back(GAP);
    pos_idx.push_back(5);
    pos_idx.push_back(6);
    pos_idx.push_back(GAP);
    pos_idx.push_back(GAP);
    pos_idx.push_back(7);
    exporter.export_pos_idx(pos_idx, oss);
    s = "6-8";
    EXPECT_EQ(s, oss.str());
}
