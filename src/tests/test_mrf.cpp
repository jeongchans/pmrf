#include <gtest/gtest.h>

#include "mrf.h"

using std::vector;
using std::pair;
using std::string;

class MRF_NodeElement_Test : public testing::Test {
  protected:
    AminoAcid abc;
};

TEST_F(MRF_NodeElement_Test, test_operator_eq) {
    MRF::NodeElement node(10);
    node.set_weight(0);
    ASSERT_EQ(10, node.get_weight().size());
    node = MRF::NodeElement(20);
    node.set_weight(0);
    ASSERT_EQ(20, node.get_weight().size());
}

TEST_F(MRF_NodeElement_Test, test_set_and_get_weight) {
    size_t num_var = abc.get_canonical_size();

    MRF::NodeElement node(abc.get_canonical_size());
    node.set_weight(0);
    Float1dArray w = node.get_weight().copy();
    ASSERT_EQ(num_var, (size_t)w.size());
    EXPECT_TRUE(all(0 == node.get_weight()));

    w(3) = 0.1;
    EXPECT_TRUE(all(0 == node.get_weight()));
    node.set_weight(w);
    EXPECT_TRUE(all(w == node.get_weight()));
    w = 0;
    EXPECT_FALSE(all(w == node.get_weight()));
}

class MRF_EdgeElement_Test : public testing::Test {
  protected:
    AminoAcid abc;
};

TEST_F(MRF_EdgeElement_Test, test_operator_eq) {
    MRF::EdgeElement edge(10, 10);
    edge.set_weight(0);
    ASSERT_EQ(10, edge.get_weight().rows());
    ASSERT_EQ(10, edge.get_weight().cols());
    edge = MRF::EdgeElement(20, 20);
    edge.set_weight(0);
    ASSERT_EQ(20, edge.get_weight().rows());
    ASSERT_EQ(20, edge.get_weight().cols());
}

TEST_F(MRF_EdgeElement_Test, test_set_and_get_weight) {
    size_t num_var = abc.get_canonical_size();

    MRF::EdgeElement edge(abc.get_canonical_size(), abc.get_canonical_size());
    edge.set_weight(0);
    Float2dArray w = edge.get_weight().copy();
    ASSERT_EQ(num_var, (size_t)w.rows());
    ASSERT_EQ(num_var, (size_t)w.cols());
    EXPECT_TRUE(all(0 == edge.get_weight()));

    w(2, 3) = 0.1;
    EXPECT_TRUE(all(0 == edge.get_weight()));
    edge.set_weight(w);
    EXPECT_TRUE(all(w == edge.get_weight()));
    w = 0;
    EXPECT_FALSE(all(w == edge.get_weight()));
}

class MRF_Test : public testing::Test {
  protected:
    virtual void SetUp() {
        length = 4;
        seq = "LMWQ";
        eidxs.push_back(EdgeIndex(0, 1));
        eidxs.push_back(EdgeIndex(0, 2));
        eidxs.push_back(EdgeIndex(1, 3));
    }

    AminoAcid abc;
    size_t length;
    string seq;
    EdgeIndexVector eidxs;
};

TEST_F(MRF_Test, test_construct_by_length) {
    MRF mrf = MRF(length, abc);
    EXPECT_EQ(length, mrf.get_length());
    EXPECT_EQ("XXXX", mrf.get_seq());
    EXPECT_EQ(length * (length - 1) / 2, mrf.get_edge_idxs().size());
}

TEST_F(MRF_Test, test_construct_by_seq) {
    MRF mrf = MRF(seq, abc);
    EXPECT_EQ(length, mrf.get_length());
    EXPECT_EQ(seq, mrf.get_seq());
    EXPECT_EQ(length * (length - 1) / 2, mrf.get_edge_idxs().size());
}

TEST_F(MRF_Test, test_construct_by_length_with_edge_idxs) {
    MRF mrf = MRF(length, abc, &eidxs);
    EdgeIndexVector ret = mrf.get_edge_idxs();
    EXPECT_EQ(eidxs.size(), ret.size());
    for (size_t i = 0; i < eidxs.size(); ++i)
        EXPECT_TRUE(eidxs[i] == ret[i]);
}

TEST_F(MRF_Test, test_construct_by_seq_with_edge_idxs) {
    MRF mrf = MRF(seq, abc, &eidxs);
    EdgeIndexVector ret = mrf.get_edge_idxs();
    EXPECT_EQ(eidxs.size(), ret.size());
    for (size_t i = 0; i < eidxs.size(); ++i)
        EXPECT_TRUE(eidxs[i] == ret[i]);
}

TEST_F(MRF_Test, test_get_var_symbol) {
    MRF mrf = MRF(seq, abc);
    EXPECT_EQ(abc.get_canonical(), mrf.get_var_symbol());
}

TEST_F(MRF_Test, test_get_num_var) {
    MRF mrf = MRF(seq, abc);
    EXPECT_EQ(abc.get_canonical_size(), mrf.get_num_var());
}

TEST_F(MRF_Test, test_get_node) {
    MRF mrf = MRF(length, abc);
    mrf.get_node(0).set_weight(0);
    Float1dArray w = mrf.get_node(0).get_weight().copy();
    ASSERT_EQ(abc.get_canonical_size(), (size_t)w.size());
    EXPECT_TRUE(all(0 == w));

    w(3) = 0.1;
    EXPECT_FALSE(all(w == mrf.get_node(0).get_weight()));
    mrf.get_node(0).set_weight(w);
    EXPECT_TRUE(all(w == mrf.get_node(0).get_weight()));
}

TEST_F(MRF_Test, test_get_edge) {
    MRF mrf = MRF(length, abc);
    mrf.get_edge(0, 2).set_weight(0);
    Float2dArray w = mrf.get_edge(0, 2).get_weight().copy();
    ASSERT_EQ(abc.get_canonical_size(), (size_t)w.rows());
    ASSERT_EQ(abc.get_canonical_size(), (size_t)w.cols());
    EXPECT_TRUE(all(0 == w));

    w(2, 3) = 0.1;
    EXPECT_FALSE(all(w == mrf.get_edge(0, 2).get_weight()));
    mrf.get_edge(0, 2).set_weight(w);
    EXPECT_TRUE(all(w == mrf.get_edge(0, 2).get_weight()));
}

TEST_F(MRF_Test, test_traverse) {
    // Custom visitor for testing
    class MyVisitor : public MRF::Visitor {
      public:
        virtual void visit_node(MRF::NodeElement*, const size_t& idx)
            { nvisits.push_back(idx); }
        virtual void visit_edge(MRF::EdgeElement*, const size_t& idx1, const size_t& idx2)
            { evisits.push_back(pair<size_t, size_t>(idx1, idx2)); }

        vector<size_t> nvisits;
        vector<pair<size_t, size_t> > evisits;
    };

    MyVisitor visitor;
    MRF mrf = MRF(length, abc);
    mrf.traverse(visitor);

    vector<size_t>::iterator npos = visitor.nvisits.begin();
    EXPECT_TRUE(0 == *npos); ++npos;
    EXPECT_TRUE(1 == *npos); ++npos;
    EXPECT_TRUE(2 == *npos); ++npos;
    EXPECT_TRUE(3 == *npos); ++npos;
    EXPECT_TRUE(npos == visitor.nvisits.end());

    vector<pair<size_t, size_t> >::iterator epos = visitor.evisits.begin();
    EXPECT_TRUE(0 == epos->first && 1 == epos->second); ++epos;
    EXPECT_TRUE(0 == epos->first && 2 == epos->second); ++epos;
    EXPECT_TRUE(0 == epos->first && 3 == epos->second); ++epos;
    EXPECT_TRUE(1 == epos->first && 2 == epos->second); ++epos;
    EXPECT_TRUE(1 == epos->first && 3 == epos->second); ++epos;
    EXPECT_TRUE(2 == epos->first && 3 == epos->second); ++epos;
    EXPECT_TRUE(epos == visitor.evisits.end());
}
