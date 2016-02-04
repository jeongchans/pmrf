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

class MRFExporter_Test : public testing::Test {
  protected:
    virtual void SetUp() {
        s = "";
        oss.flush();
    }

    MRFExporter exporter;
    string s;
    std::ostringstream oss;
};

TEST_F(MRFExporter_Test, test_export_seq) {
    string seq = "PPDQEFLRARVQLGDA";
    s = "PPDQE\nFLRAR\nVQLGD\nA\n";
    exporter.export_seq(seq, 5, oss);
    EXPECT_EQ(s, oss.str());
}

TEST_F(MRFExporter_Test, test_export_node_symbol) {
    string sym = "LWQ";
    exporter.export_node_symbol(sym, oss);
    s = "L\tW\tQ";
    EXPECT_EQ(s, oss.str());
}

TEST_F(MRFExporter_Test, test_export_node_weight) {
    Float1dArray w(3);
    w = 0.1, 0.0, 0.9;
    s = "0.1\t*\t0.9";
    exporter.export_node_weight(w, oss);
    EXPECT_EQ(s, oss.str());
}

TEST_F(MRFExporter_Test, test_export_edge_symbol) {
    string sym = "LWQ";
    exporter.export_edge_symbol(sym, oss);
    s = "LL\tLW\tLQ\tWL\tWW\tWQ\tQL\tQW\tQQ";
    EXPECT_EQ(s, oss.str());
}

TEST_F(MRFExporter_Test, test_export_edge_weight) {
    Float2dArray w(3, 3);
    w = 0.0;
    w(0, 1) = 0.1;
    w(1, 0) = 0.5;
    w(2, 2) = 0.4;
    s = "*\t0.1\t*\t0.5\t*\t*\t*\t*\t0.4";
    exporter.export_edge_weight(w, oss);
    EXPECT_EQ(s, oss.str());
}
