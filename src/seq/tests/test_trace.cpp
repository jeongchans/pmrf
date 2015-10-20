#include <gtest/gtest.h>

#include "trace.h"

// Example: PPDQEFLRG---ARVQLGDA
//          --DQ---HGNRIVHLQ----
string seq1  = "PPDQEFLRGARVQLGDA";
string seq2  = "DQHGNRIVHLQ";
string aseq1 = "PPDQEFLRGARVQLGDA";
string aseq2 = "--DQ---HGNRIVHLQ----";
string st1   = "MMMMMMMMMMMMMMMMM";
//string st2   = "UUMMOEEMMIIIMMMMUUUU";
string st2   = "DDMMDDDMMIIIMMMMDDDD";

class TraceTest : public testing::Test {
  protected:
    virtual void SetUp() {
        trace1 = Trace(st1, seq1);
        trace2 = Trace(st2, seq2);
    }

    Trace trace1, trace2;
};

TEST_F(TraceTest, test_get_seq) {
    EXPECT_EQ(seq1, trace1.get_seq());
    EXPECT_EQ(seq2, trace2.get_seq());
}

TEST_F(TraceTest, test_is_matched) {
    EXPECT_TRUE(trace1.is_matched(3));
    EXPECT_TRUE(trace2.is_matched(3));
    EXPECT_FALSE(trace2.is_matched(4));
    EXPECT_FALSE(trace2.is_matched(6));
    EXPECT_TRUE(trace2.is_matched(7));
}

TEST_F(TraceTest, test_get_trimmed_seq) {
    EXPECT_EQ(string("PPDQEFLRGARVQLGDA"), trace1.get_trimmed_aseq());
    //EXPECT_EQ(string("^^DQ=--HGVHLQ^^^^"), trace2.get_trimmed_aseq());
    EXPECT_EQ(string("--DQ---HGVHLQ----"), trace2.get_trimmed_aseq());
}

TEST_F(TraceTest, test_operator_eq) {
    Trace trace(st1, seq1);
    EXPECT_TRUE(trace1 == trace);
    EXPECT_FALSE(trace2 == trace);
}

TEST_F(TraceTest, test_get_symbol_at) {
    EXPECT_EQ('P', trace1.get_symbol_at(0));
    //EXPECT_EQ('^', trace2.get_symbol_at(0));
    EXPECT_EQ('-', trace2.get_symbol_at(0));
    EXPECT_EQ('H', trace2.get_symbol_at(7));
}

class TraceVectorTest : public testing::Test {
  public:
    TraceVectorTest() {
        traces.push_back(Trace(st1, seq1));
        traces.push_back(Trace(st2, seq2));
    }

  protected:
    AminoAcid abc;
    TraceVector traces;
};

TEST_F(TraceVectorTest, test_subset_matched) {
    TraceVector trs = traces.subset_matched(0);
    ASSERT_EQ((size_t) 1, trs.size());
    EXPECT_TRUE(traces[0] == trs[0]);

    trs = traces.subset_matched(2);
    ASSERT_EQ((size_t) 2, trs.size());
}

TEST_F(TraceVectorTest, test_get_trimmed_aseq_vec) {
    vector<string> vec = traces.get_trimmed_aseq_vec();
    ASSERT_EQ(traces.size(), vec.size());
    EXPECT_EQ(traces[0].get_trimmed_aseq(), vec[0]);
    EXPECT_EQ(traces[1].get_trimmed_aseq(), vec[1]);
}

string a3m(">seq1 Sample #1\n"
           "PPDQEFLRGARVQLGDA\n"
           ">seq2 Sample #2\n"
           "--DQ---HGnriVHLQ----\n"
           ">seq3 Ignored\n"
           "-----------------agv\n");

string afa(">seq1 Sample #1\n"
           "PPDQEFLRG---ARVQLGDA---\n"
           ">seq2 Sample #2\n"
           "--DQ---HGNRIVHLQ-------\n"
           ">seq3 Ignored\n"
           "--------------------AGV\n");

class TraceImporterTest : public testing::Test {
  public:
    TraceImporterTest() {
        traces.push_back(Trace(st1, seq1));
        traces.push_back(Trace(st2, seq2));
    }

  protected:
    AminoAcid abc;
    TraceVector traces;
};

TEST_F(TraceImporterTest, test_import_a3m) {
    std::istringstream is(a3m);
    TraceImporter importer(abc);
    std::pair<size_t, TraceVector> ret = importer.import(is, A3M);
    size_t& length = ret.first;
    TraceVector& trs = ret.second;
    ASSERT_EQ((size_t) 17, length);
    ASSERT_EQ(traces.size(), trs.size());
    for (size_t i = 0; i < traces.size(); ++i) {
        EXPECT_TRUE(traces[i] == trs[i]) 
        << "i = " << i << endl                                    
        << trs[i].get_seq() << endl
        << trs[i].get_trimmed_aseq() << endl;
    }
}

TEST_F(TraceImporterTest, test_import_afa) {
    std::istringstream is(afa);
    TraceImporter importer(abc);
    std::pair<size_t, TraceVector> ret = importer.import(is, AFASTA);
    size_t& length = ret.first;
    TraceVector& trs = ret.second;
    ASSERT_EQ(traces.size(), trs.size());
    ASSERT_EQ((size_t) 17, length);
    for (size_t i = 0; i < traces.size(); ++i) {
        EXPECT_TRUE(traces[i] == trs[i]);
    }
}
