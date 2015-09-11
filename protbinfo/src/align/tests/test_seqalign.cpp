#include <gtest/gtest.h>

#include "seqalign.h"

class SeqAlignTest : public testing::Test {
  public:
    SeqAlignTest() {
        seq1 = "NGPSTKDFGKISESREFDNQNGPSTKDFGKISESREFDNQ";
        seq2 = "QNQLERSFGKINMRLEDALVQNQLERSFGKINMRLEDALV";
    }

    string seq1;
    string seq2;
    double gap_open;
    double gap_ext;
    double baseline;
    BLOSUM62Matrix blosum62;
};

TEST_F(SeqAlignTest, test_global_align) {
    gap_open = 7.0;
    gap_ext = 0.5;
    baseline = 0.0;
    NWAlign align(blosum62, gap_open, gap_ext, baseline);
    PairwiseAlignment aln = align(seq1, seq2);
    EXPECT_EQ("NGPSTKDFGKIS---ESREFDNQNGPSTKDFGKISESREFDNQ", aln.alignment[0].seq);
    EXPECT_EQ(0, aln.alignment[0].pos_idx.front());
    EXPECT_EQ((int)seq1.size() - 1, aln.alignment[0].pos_idx.back());
    EXPECT_EQ("QNQLERSFGKINMRLEDALVQNQ---LERSFGKINMRLEDALV", aln.alignment[1].seq);
    EXPECT_EQ(0, aln.alignment[1].pos_idx.front());
    EXPECT_EQ((int)seq2.size() - 1, aln.alignment[1].pos_idx.back());
    EXPECT_EQ(26, aln.property["score"]);
}

TEST_F(SeqAlignTest, test_local_align) {
    gap_open = 5.0;
    gap_ext = 0.5;
    baseline = 0.0;
    SWAlign align(blosum62, gap_open, gap_ext, baseline);
    PairwiseAlignment aln = align(seq1, seq2);
    EXPECT_EQ("KDFGKISESREFD----NQNGPSTKDFGKISESREFD", aln.alignment[0].seq);
    EXPECT_EQ(5, aln.alignment[0].pos_idx.front());
    EXPECT_EQ(37, aln.alignment[0].pos_idx.back());
    EXPECT_EQ("RSFGKINMRLE-DALVQNQ---LERSFGKINMRLE-D", aln.alignment[1].seq);
    EXPECT_EQ(5, aln.alignment[1].pos_idx.front());
    EXPECT_EQ(36, aln.alignment[1].pos_idx.back());
    EXPECT_EQ(45.5, aln.property["score"]);
}
