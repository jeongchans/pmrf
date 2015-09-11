#include <gtest/gtest.h>

#include "hmmio.h"

class HHsearchHMMExporterTest : public testing::Test {
  protected:
    virtual void SetUp() {
        s = "";
        oss.flush();
    }

    HHsearchHMMExporter exporter;
    std::string s;
    std::ostringstream oss;
};

TEST_F(HHsearchHMMExporterTest, test_export_emit_symbol) {
    std::string sym = "LWQP";
    exporter.export_emit_symbol(sym, oss);
    s = "L\tW\tQ\tP\t";
    EXPECT_EQ(s, oss.str());
}

TEST_F(HHsearchHMMExporterTest, test_export_emit) {
    Float1dArray em(4);
    em = 0.1, 0.5, 0.0, 0.4;
    s = "3322\t1000\t*\t1322\t";
    exporter.export_emit(em, oss);
    EXPECT_EQ(s, oss.str());
}

TEST_F(HHsearchHMMExporterTest, test_export_transit) {
    Float1dArray tr(3);
    tr = 0.7, 0.2, 0.1;
    exporter.export_transit(MATCH, tr, oss);
    tr = 0.7, 0.0, 0.3;
    exporter.export_transit(INSERT, tr, oss);
    tr = 0.6, 0.4, 0.0;
    exporter.export_transit(DELETE, tr, oss);
    s = "515\t3322\t2322\t515\t1737\t737\t1322\t";
    EXPECT_EQ(s, oss.str());
}

TEST_F(HHsearchHMMExporterTest, test_export_eff_num) {
    double neff;
    neff = 4.5;
    exporter.export_eff_num(neff, oss);
    neff = 1;
    exporter.export_eff_num(neff, oss);
    neff = 0;
    exporter.export_eff_num(neff, oss);
    s = "2459\t1000\t0\t";
    EXPECT_EQ(s, oss.str());
}

class HHsearchHMMImporterTest : public testing::Test {
  public:
    virtual void SetUp() {
        s =  "";
    }

    std::string s;
    AminoAcid abc;
};

TEST_F(HHsearchHMMImporterTest, test_import_emit_symbol) {
    HHsearchHMMImporter importer(abc);
    s = "L\tW\tQ\tP\t";
    std::istringstream is(s);
    EXPECT_EQ("LWQP", importer.import_emit_symbol(is));
}

TEST_F(HHsearchHMMImporterTest, test_import_emit) {
    HHsearchHMMImporter importer(abc);
    s = "3322\t1000\t*\t1322\t";
    std::istringstream is(s);
    int size = 4;
    Float1dArray em(size);
    em = 0.1, 0.5, 0.0, 0.4;
    EXPECT_TRUE(allclose(em, importer.import_emit(is, size)));
}

TEST_F(HHsearchHMMImporterTest, test_import_transit) {
    HHsearchHMMImporter importer(abc);
    s = "515\t3322\t2322\t515\t1737\t737\t1322\t";
    std::istringstream is(s);
    Float1dArray tr(3);
    tr = 0.7, 0.2, 0.1;
    EXPECT_TRUE(allclose(tr, importer.import_transit(is, MATCH)));
    tr = 0.7, 0.0, 0.3;
    EXPECT_TRUE(allclose(tr, importer.import_transit(is, INSERT)));
    tr = 0.6, 0.4, 0.0;
    EXPECT_TRUE(allclose(tr, importer.import_transit(is, DELETE)));
}

TEST_F(HHsearchHMMImporterTest, test_import_eff_num) {
    HHsearchHMMImporter importer(abc);
    s = "2459\t1000\t0\t";
    std::istringstream is(s);
    double n = importer.import_eff_num(is);
    EXPECT_TRUE(allclose(4.5, n, 1e-2)) << n;
    EXPECT_TRUE(allclose(1, importer.import_eff_num(is)));
    EXPECT_TRUE(allclose(0, importer.import_eff_num(is)));
}

TEST_F(HHsearchHMMImporterTest, test_consistency) {
    ProfileHMM model(3, abc);
    TraceVector traces;
    HMMParamEstimator estimator(traces, abc);
    estimator.parameterize(model);
    Float1dArray tr(3);
    tr = 0.7, 0.2, 0.1;
    model.get_match(1).set_transit(tr);

    std::ostringstream oss;
    HHsearchHMMExporter exporter;
    exporter.export_model(model, oss);
    s = oss.str();

    std::istringstream iss(s);
    HHsearchHMMImporter importer(abc);
    ProfileHMM reload = importer.import_model(iss);
    EXPECT_TRUE(allclose(tr, reload.get_match(1).get_transit()));

    std::ostringstream oss2;
    exporter.export_model(reload, oss2);
    EXPECT_EQ(s, oss2.str());
}

//TEST_F(HCIModelImporterTest, test_extract_hcisvm_content) {
//    HCIModelImporter importer(abc);
//    string s = "line 1\nline2\nline3\nHCISVM_END\nline4";
//    std::istringstream is(s);
//    string content = importer.extract_hcisvm_content(is);
//    EXPECT_EQ("line 1\nline2\nline3\n", content) << content;
//}
