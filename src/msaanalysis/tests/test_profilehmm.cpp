#include <gtest/gtest.h>

#include "profilehmm.h"

class ProfileHMMTest : public testing::Test {
  protected:
    virtual void SetUp() {
        length = 4;
        seq = "LMWQ";
    }

    AminoAcid abc;
    size_t length;
    std::string seq;
};

TEST_F(ProfileHMMTest, test_construct_by_length) {
    ProfileHMM hmm = ProfileHMM(length, abc);
    EXPECT_EQ(length, hmm.get_length());
    EXPECT_EQ("XXXX", hmm.get_seq());
}

TEST_F(ProfileHMMTest, test_construct_by_seq) {
    ProfileHMM hmm = ProfileHMM(seq, abc);
    EXPECT_EQ(length, hmm.get_length());
    EXPECT_EQ(seq, hmm.get_seq());
}

TEST_F(ProfileHMMTest, test_get_emit_symbol) {
    ProfileHMM hmm = ProfileHMM(length, abc);
    EXPECT_EQ(abc.get_canonical(), hmm.get_emit_symbol());
}

TEST_F(ProfileHMMTest, test_get_emit) {
    ProfileHMM hmm = ProfileHMM(length, abc);
    Float1dArray v = hmm.get_match(0).get_emit().copy();
    ASSERT_EQ(abc.get_canonical_size(), (size_t)v.size());
    EXPECT_TRUE(all(0 == v));

    v = 0.1, 0.6, 0.2, 0.1;
    EXPECT_FALSE(all(v == hmm.get_match(0).get_emit()));
    hmm.get_match(0).set_emit(v);
    EXPECT_TRUE(all(v == hmm.get_match(0).get_emit()));
}

TEST_F(ProfileHMMTest, test_get_transit) {
    ProfileHMM hmm = ProfileHMM(length, abc);
    Float1dArray v = hmm.get_match(0).get_transit().copy();
    ASSERT_EQ(NUM_STATE_TYPE, (size_t)v.size());
    EXPECT_TRUE(all(0 == v));

    v = 0.3, 0.4, 0.3;
    EXPECT_FALSE(all(v == hmm.get_match(0).get_transit()));
    hmm.get_match(0).set_transit(v);
    EXPECT_TRUE(all(v == hmm.get_match(0).get_transit()));
}

TEST_F(ProfileHMMTest, test_traverse_visitor) {
    // Custom visitor
    typedef std::pair<StateType, size_t> Node;
    class MyVisitor : public HMMStateVisitor {
      public:
        virtual void visit_match(HMMMatchState*, const size_t& idx)
            { add_visit(MATCH, idx); }
        virtual void visit_delete(HMMDeleteState*, const size_t& idx)
            { add_visit(DELETE, idx); }
        virtual void visit_insert(HMMInsertState*, const size_t& idx)
            { add_visit(INSERT, idx); }

        std::vector<Node> visits;

      private:
        void add_visit(const StateType& type, const size_t& idx) {
            visits.push_back(Node(type, idx));
        }
    };

    MyVisitor visitor;
    ProfileHMM hmm = ProfileHMM(length, abc);
    hmm.traverse(visitor);
    size_t tot_size = length * NUM_STATE_TYPE;
    ASSERT_EQ(tot_size, visitor.visits.size());
    std::vector<Node>::iterator pos = visitor.visits.begin();
    EXPECT_TRUE(MATCH == pos->first  && 0 == pos->second);  ++pos;
    EXPECT_TRUE(DELETE == pos->first && 0 == pos->second);  ++pos;
    EXPECT_TRUE(INSERT == pos->first && 0 == pos->second);  ++pos;
    EXPECT_TRUE(MATCH == pos->first  && 1 == pos->second);  ++pos;
    EXPECT_TRUE(DELETE == pos->first && 1 == pos->second);  ++pos;
    EXPECT_TRUE(INSERT == pos->first && 1 == pos->second);  ++pos;
    EXPECT_TRUE(MATCH == pos->first  && 2 == pos->second);  ++pos;
    EXPECT_TRUE(DELETE == pos->first && 2 == pos->second);  ++pos;
    EXPECT_TRUE(INSERT == pos->first && 2 == pos->second);  ++pos;
    EXPECT_TRUE(MATCH == pos->first  && 3 == pos->second);  ++pos;
    EXPECT_TRUE(DELETE == pos->first && 3 == pos->second);  ++pos;
    EXPECT_TRUE(INSERT == pos->first && 3 == pos->second);  ++pos;
    EXPECT_TRUE(pos == visitor.visits.end());
}

class HMMParamEstimatorTest : public testing::Test {
  protected:
    virtual void SetUp() {
        traces.push_back(Trace("MMMMMMMMMMMMMMMM", "PPDQEFLRARVQLGDA"));
        traces.push_back(Trace("DDDDDMMMIIIIMMMMDDDD", "SSHGNRIVHLQ"));
    }

    AminoAcid abc;
    TraceVector traces;
};

TEST_F(HMMParamEstimatorTest, test_visit_match) {
    HMMParamEstimator estimator(traces, abc);
    Float1dArray em(abc.get_canonical_size());
    Float1dArray tr(NUM_STATE_TYPE);
    HMMMatchState state(abc.get_canonical_size());

    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 1.0, 0.0, 0.0;
    estimator.visit_match(&state, 0);
    EXPECT_TRUE(all(em == state.get_emit())) << state.get_emit() << std::endl;
    EXPECT_TRUE(all(tr == state.get_transit())) << state.get_transit() << std::endl;
    EXPECT_EQ(1, state.get_eff_num());

    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 0.5, 0.0, 0.5;
    estimator.visit_match(&state, 7);
    EXPECT_TRUE(all(em == state.get_emit()));
    EXPECT_TRUE(all(tr == state.get_transit()));
    EXPECT_EQ(2, state.get_eff_num());
}

TEST_F(HMMParamEstimatorTest, test_visit_insert) {
    HMMParamEstimator estimator(traces, abc);
    Float1dArray em(abc.get_canonical_size());
    Float1dArray tr(NUM_STATE_TYPE);
    HMMInsertState state(abc.get_canonical_size());

    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 0.0, 0.0, 0.0;
    estimator.visit_insert(&state, 0);
    EXPECT_TRUE(all(em == state.get_emit())) << state.get_emit() << std::endl;
    EXPECT_TRUE(all(tr == state.get_transit())) << state.get_transit() << std::endl;

    em = 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.00, 0.25, 0.00, 0.00,
         0.00, 0.25, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00;
    tr = 0.25, 0.00, 0.75;
    estimator.visit_insert(&state, 7);
    EXPECT_TRUE(all(em == state.get_emit()));
    EXPECT_TRUE(all(tr == state.get_transit()));
}

TEST_F(HMMParamEstimatorTest, test_visit_delete) {
    HMMParamEstimator estimator(traces, abc);
    Float1dArray tr(NUM_STATE_TYPE);
    HMMDeleteState state;

    tr = 0.0, 1.0, 0.0;
    estimator.visit_delete(&state, 0);
    EXPECT_TRUE(all(tr == state.get_transit())) << state.get_transit() << std::endl;

    tr = 1.0, 0.0, 0.0;
    estimator.visit_delete(&state, 4);
    EXPECT_TRUE(all(tr == state.get_transit()));
}

TEST_F(HMMParamEstimatorTest, test_termi_gap_filter) {
    HMMParamEstimator estimator(traces, abc, true);
    Float1dArray tr(NUM_STATE_TYPE);
    HMMDeleteState state;

    tr = 0.0, 0.0, 0.0;
    estimator.visit_delete(&state, 4);
    EXPECT_TRUE(all(tr == state.get_transit())) << state.get_transit() << std::endl;

    tr = 0.0, 0.0, 0.0;
    estimator.visit_delete(&state, 12);
    EXPECT_TRUE(all(tr == state.get_transit()));
}

class ParameterizationTest : public testing::Test {
  protected:
    virtual void SetUp() {
        em.resize(abc.get_canonical_size());
        tr.resize(NUM_STATE_TYPE);
        length = 16;
        // training data
        traces.push_back(Trace("MMMMMMMMMMMMMMMM", "PPDQEFLRARVQLGDA"));
        traces.push_back(Trace("DDDDDMMMIIIIMMMMDDDD", "SSHGNRIVHLQ"));
    }

    AminoAcid abc;
    size_t length;
    TraceVector traces;
    Float1dArray em, tr;
};

TEST_F(ParameterizationTest, test_state_accept) {
    HMMParamEstimator estimator(traces, abc);  // estimator setting

    HMMMatchState m_state(abc.get_canonical_size());
    m_state.accept(&estimator, 0);
    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 1.0, 0.0, 0.0;
    EXPECT_TRUE(all(em == m_state.get_emit())) << m_state.get_emit() << std::endl;
    EXPECT_TRUE(all(tr == m_state.get_transit())) << m_state.get_transit() << std::endl;

    HMMInsertState i_state(abc.get_canonical_size());
    i_state.accept(&estimator, 0);
    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 0.0, 0.0, 0.0;
    EXPECT_TRUE(all(em == i_state.get_emit())) << i_state.get_emit() << std::endl;
    EXPECT_TRUE(all(tr == i_state.get_transit())) << i_state.get_transit() << std::endl;

    HMMDeleteState d_state;
    d_state.accept(&estimator, 0);
    tr = 0.0, 1.0, 0.0;
    EXPECT_TRUE(all(tr == d_state.get_transit())) << d_state.get_transit() << std::endl;
}

TEST_F(ParameterizationTest, test_no_strategy) {
    ProfileHMM hmm = ProfileHMM(length, abc);
    HMMParamEstimator estimator(traces, abc);  // estimator setting
    estimator.parameterize(hmm);                        // parameterization

    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 1.0, 0.0, 0.0;
    EXPECT_TRUE(all(em == hmm.get_match(0).get_emit())) << hmm.get_match(0).get_emit() << std::endl;
    EXPECT_TRUE(all(tr == hmm.get_match(0).get_transit())) << hmm.get_match(0).get_transit() << std::endl;

    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 0.5, 0.0, 0.5;
    EXPECT_TRUE(all(em == hmm.get_match(7).get_emit())) << hmm.get_match(7).get_emit() << std::endl;
    EXPECT_TRUE(all(tr == hmm.get_match(7).get_transit())) << hmm.get_match(7).get_transit() << std::endl;

    em = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    tr = 0.0, 0.0, 0.0;
    EXPECT_TRUE(all(em == hmm.get_insert(0).get_emit()));
    EXPECT_TRUE(all(tr == hmm.get_insert(0).get_transit()));

    em = 0.00, 0.00, 0.00, 0.00, 0.00, 0.25, 0.00, 0.25, 0.00, 0.00,
         0.00, 0.25, 0.00, 0.00, 0.25, 0.00, 0.00, 0.00, 0.00, 0.00;
    tr = 0.25, 0.00, 0.75;
    EXPECT_TRUE(all(em == hmm.get_insert(7).get_emit()));
    EXPECT_TRUE(all(tr == hmm.get_insert(7).get_transit()));

    tr = 0.0, 1.0, 0.0;
    EXPECT_TRUE(all(tr == hmm.get_delete(0).get_transit()));

    tr = 1.0, 0.0, 0.0;
    EXPECT_TRUE(all(tr == hmm.get_delete(4).get_transit()));
}
