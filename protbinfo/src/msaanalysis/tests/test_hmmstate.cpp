#include <gtest/gtest.h>

#include "hmmstate.h"

class HMMStateTest : public testing::Test {
  protected:
    virtual void SetUp() {
        m_state = HMMMatchState(abc.get_canonical_size());
        d_state = HMMDeleteState();
        i_state = HMMInsertState(abc.get_canonical_size());
    }

    AminoAcid abc;
    HMMMatchState m_state;
    HMMDeleteState d_state;
    HMMInsertState i_state;
};

TEST_F(HMMStateTest, test_set_and_get_transit) {
    Float1dArray v = m_state.get_transit().copy();
    ASSERT_EQ(NUM_STATE_TYPE, (size_t)v.size());
    EXPECT_TRUE(all(0 == m_state.get_transit()));
    EXPECT_TRUE(all(0 == d_state.get_transit()));
    EXPECT_TRUE(all(0 == i_state.get_transit()));

    v = 0.3, 0.3, 0.4;
    EXPECT_TRUE(all(0 == m_state.get_transit()));
    m_state.set_transit(v);
    EXPECT_TRUE(all(v == m_state.get_transit()));
    EXPECT_TRUE(all(0 == d_state.get_transit()));
    EXPECT_TRUE(all(0 == i_state.get_transit()));

    EXPECT_FLOAT_EQ(0.3, m_state.get_transit_to(MATCH));
    EXPECT_FLOAT_EQ(0.3, m_state.get_transit_to(DELETE));
    EXPECT_FLOAT_EQ(0.4, m_state.get_transit_to(INSERT));
}

TEST_F(HMMStateTest, test_set_and_get_eff_num) {
    EXPECT_EQ(0, m_state.get_eff_num());

    double neff = 2.3;
    m_state.set_eff_num(neff);
    EXPECT_EQ(neff, m_state.get_eff_num());
}

TEST_F(HMMStateTest, test_operator_eq) {
    HMMMatchState state(3);
    ASSERT_EQ(3, state.get_emit().size());
    state = HMMMatchState(5);
    ASSERT_EQ(5, state.get_emit().size());
}

TEST_F(HMMStateTest, test_set_and_get_emit) {
    Float1dArray v = m_state.get_emit().copy();
    ASSERT_EQ(abc.get_canonical_size(), (size_t)v.size());
    EXPECT_TRUE(all(0 == m_state.get_emit()));
    EXPECT_TRUE(all(0 == i_state.get_emit()));

    v(3) = 0.1;
    EXPECT_TRUE(all(0 == m_state.get_emit()));
    m_state.set_emit(v);
    EXPECT_TRUE(all(v == m_state.get_emit()));
    EXPECT_TRUE(all(0 == i_state.get_emit()));
    v = 0;
    EXPECT_FALSE(all(v == m_state.get_emit()));
}

TEST_F(HMMStateTest, test_set_and_get_property) {
    Float1dArray v(3);
    v = 1, 0, 2;
    m_state.set_property("v", v);
    Float1dArray r_v = *(Float1dArray*)(m_state.get_property("v"));
    EXPECT_TRUE(all(v == r_v));

    FloatType f = 0.5;
    m_state.set_property("f", f);
    FloatType r_f = *(FloatType*)(m_state.get_property("f"));
    EXPECT_FLOAT_EQ(f, r_f);

    Float2dArray m(2, 3);
    m = 0, 1, 2,
        3, 4, 5;
    m_state.set_property("m", m);
    Float2dArray r_m = *(Float2dArray*)(m_state.get_property("m"));
    EXPECT_TRUE(all(m == r_m));
}
