#include <gtest/gtest.h>

#include "trace.h"

class TraceTest : public testing::Test {
  protected:
    virtual void SetUp() {
        st1 = "MMMMMMMMMMMMMMMM";
        seq1 = "PPDQEFLRARVQLGDA";
        trace1 = Trace(st1, seq1);
        st2 = "DDDDDMMMIIIIMMMMDDDD";
        seq2 = "SSHGNRIVHLQ";
        trace2 = Trace(st2, seq2);
    }

    Trace trace1, trace2;
    std::string st1, seq1, st2, seq2;
};

TEST_F(TraceTest, test_get_seq) {
    EXPECT_EQ(seq1, trace1.get_seq());
    EXPECT_EQ(seq2, trace2.get_seq());
}

TEST_F(TraceTest, test_is_passing) {
    EXPECT_TRUE(trace1.is_passing(MATCH, 3));
    EXPECT_FALSE(trace1.is_passing(INSERT, 3));
    EXPECT_FALSE(trace1.is_passing(DELETE, 3));
    EXPECT_TRUE(trace2.is_passing(DELETE, 3));
    EXPECT_FALSE(trace2.is_passing(MATCH, 3));
    EXPECT_FALSE(trace2.is_passing(INSERT, 3));
    EXPECT_TRUE(trace2.is_passing(MATCH, 7));
    EXPECT_TRUE(trace2.is_passing(INSERT, 7));
    EXPECT_TRUE(trace2.is_passing(MATCH, 8));
}

TEST_F(TraceTest, test_has_terminal_gap) {
    EXPECT_FALSE(trace1.has_terminal_gap(3));
    EXPECT_TRUE(trace2.has_terminal_gap(3));
    EXPECT_TRUE(trace2.has_terminal_gap(4));
    EXPECT_FALSE(trace2.has_terminal_gap(5));
    EXPECT_FALSE(trace2.has_terminal_gap(11));
    EXPECT_TRUE(trace2.has_terminal_gap(12));
}

TEST_F(TraceTest, test_get_emit) {
    EXPECT_EQ(std::string("P"), trace1.get_emit(MATCH, 0));
    EXPECT_EQ(std::string(""), trace1.get_emit(DELETE, 0));
    EXPECT_EQ(std::string(""), trace2.get_emit(MATCH, 0));
    EXPECT_EQ(std::string("H"), trace2.get_emit(MATCH, 7));
    EXPECT_EQ(std::string("GNRI"), trace2.get_emit(INSERT, 7));
}

TEST_F(TraceTest, test_get_visit) {
    EXPECT_EQ(1, trace1.get_visit(MATCH, 3));
    EXPECT_EQ(0, trace1.get_visit(DELETE, 3));
    EXPECT_EQ(0, trace1.get_visit(INSERT, 3));
    EXPECT_EQ(0, trace2.get_visit(MATCH, 3));
    EXPECT_EQ(1, trace2.get_visit(DELETE, 3));
    EXPECT_EQ(0, trace2.get_visit(INSERT, 3));
    EXPECT_EQ(1, trace2.get_visit(MATCH, 7));
    EXPECT_EQ(4, trace2.get_visit(INSERT, 7));
}

TEST_F(TraceTest, test_get_transit_count) {
    FreqVec cnt(NUM_STATE_TYPE);
    cnt = 0, 0, 0;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(MATCH, 0)));
    cnt = 0, 1, 0;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(DELETE, 0)));
    cnt = 1, 0, 0;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(DELETE, 4)));
    cnt = 1, 0, 0;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(MATCH, 5)));
    cnt = 0, 0, 1;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(MATCH, 7)));
    cnt = 1, 0, 3;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(INSERT, 7)));
    cnt = 0, 1, 0;
    EXPECT_TRUE(all(cnt == trace2.get_transit_count(MATCH, 11)));
}

TEST_F(TraceTest, test_get_MD_seq) {
    EXPECT_EQ(std::string("PPDQEFLRARVQLGDA"), trace1.get_MD_seq());
    EXPECT_EQ(std::string("-----SSHVHLQ----"), trace2.get_MD_seq());
}

TEST_F(TraceTest, test_operator_eq) {
    Trace trace("MMMMMMMMMMMMMMMM",
                "PPDQEFLRARVQLGDA");
    EXPECT_TRUE(trace1 == trace);
    EXPECT_FALSE(trace2 == trace);
}

class TraceVectorTest : public testing::Test {
  public:
    TraceVectorTest() {
        traces.push_back(Trace("MMMMMMMMMMMMMMMM", "PPDQEFLRARVQLGDA"));
        traces.push_back(Trace("DDDDDMMMIIIIMMMMDDDD", "SSHGNRIVHLQ"));
    }

  protected:
    AminoAcid abc;
    TraceVector traces;
};

TEST_F(TraceVectorTest, test_subset_passing) {
    TraceVector trs = traces.subset_passing(MATCH, 0);
    ASSERT_EQ((size_t) 1, trs.size());
    EXPECT_TRUE(traces[0] == trs[0]);

    trs = traces.subset_passing(DELETE, 0, true);
    EXPECT_TRUE(trs.empty());
}

TEST_F(TraceVectorTest, test_get_MD_seq_vec) {
    vector<string> vec = traces.get_MD_seq_vec();
    ASSERT_EQ(traces.size(), vec.size());
    EXPECT_EQ(traces[0].get_MD_seq(), vec[0]);
    EXPECT_EQ(traces[1].get_MD_seq(), vec[1]);
}

class TraceImporterTest : public testing::Test {
  public:
    TraceImporterTest() {
        traces.push_back(Trace("MMMMMMMMMMMMMMMM", "PPDQEFLRARVQLGDA"));
        traces.push_back(Trace("DDDDDMMMIIIIMMMMDDDD", "SSHGNRIVHLQ"));
    }

  protected:
    AminoAcid abc;
    TraceVector traces;
};

TEST_F(TraceImporterTest, test_import_a3m) {
    std::string a3m(
        ">seq1 Sample #1\n"
        "PPDQEFLRARVQLGDA\n"
        ">seq2 Sample #2\n"
        "-----SSHgnriVHLQ----\n"
        ">seq3 Ignored\n"
        "----------------agv\n"
        );
    std::istringstream is(a3m);

    TraceImporter importer(abc);
    std::pair<size_t, TraceVector> ret = importer.import(is, A3M);
    size_t& length = ret.first;
    TraceVector& trs = ret.second;
    ASSERT_EQ(traces.size(), trs.size());
    ASSERT_EQ((size_t) 16, length);
    for (size_t i = 0; i < traces.size(); ++i) {
        EXPECT_TRUE(traces[i] == trs[i]);
    }
}

TEST_F(TraceImporterTest, test_import_afa) {
    std::string afa(
        ">seq1 Sample #1\n"
        "PPDQEFLR----ARVQLGDA---\n"
        ">seq2 Sample #2\n"
        "-----SSHGNRIVHLQ-------\n"
        ">seq3 Ignored\n"
        "--------------------AGV\n"
        );
    std::istringstream is(afa);

    TraceImporter importer(abc);
    std::pair<size_t, TraceVector> ret = importer.import(is, AFASTA);
    size_t& length = ret.first;
    TraceVector& trs = ret.second;
    ASSERT_EQ(traces.size(), trs.size());
    ASSERT_EQ((size_t) 16, length);
    for (size_t i = 0; i < traces.size(); ++i) {
        EXPECT_TRUE(traces[i] == trs[i]);
    }
}

TEST_F(TraceImporterTest, test_import_afa_with_gap) {
    const Alphabet aa_gap("ACDEFGHIKLMNPQRSTVWY-",
                          "",
                          "BJZOU",
                          "X",
                          "*",
                          "~");
    std::string afa(
        ">seq1 Sample #1\n"
        "PPDQEFLR----ARVQLGDA---\n"
        ">seq2 Sample #2\n"
        "-----SSHGNRIVHLQ-------\n"
        ">seq3 Ignored\n"
        "--------------------AGV\n"
        );
    std::istringstream is(afa);
    TraceVector trvec;
    trvec.push_back(Trace("MMMMMMMMMMMMMMMM", "PPDQEFLRARVQLGDA"));
    trvec.push_back(Trace("MMMMMMMMIIIIMMMMMMMM", "-----SSHGNRIVHLQ----"));

    TraceImporter importer(aa_gap);
    std::pair<size_t, TraceVector> ret = importer.import(is, AFASTA);
    size_t& length = ret.first;
    TraceVector& trs = ret.second;
    ASSERT_EQ(trvec.size(), trs.size());
    ASSERT_EQ((size_t) 16, length);
    for (size_t i = 0; i < trvec.size(); ++i) {
        EXPECT_TRUE(trvec[i] == trs[i]);
    }
}
