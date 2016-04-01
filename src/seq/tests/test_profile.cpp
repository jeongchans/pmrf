#include <gtest/gtest.h>

#include "profile.h"

class ProfileTest : public testing::Test {
  public:
    ProfileTest() : length(4), seq("LMWQ") {};

    AminoAcid abc;
    size_t length;
    std::string seq;
};

TEST_F(ProfileTest, test_construct_by_length) {
    Profile profile = Profile(length, abc);
    EXPECT_EQ(length, profile.get_length());
    EXPECT_EQ("XXXX", profile.get_seq());
}

TEST_F(ProfileTest, test_construct_by_seq) {
    Profile profile = Profile(seq, abc);
    EXPECT_EQ(length, profile.get_length());
    EXPECT_EQ(seq, profile.get_seq());
}

class ProfileBuilderTest : public testing::Test {
  public:
    ProfileBuilderTest() : length(16) {};

    virtual void SetUp() {
        traces.push_back(Trace("MMMMMMMMMMMMMMMM", "PPDQEFLRARVQLGDA"));
//        traces.push_back(Trace("UUUUUMMMIIIIMMMMUUUU", "SSHGNRIVHLQ"));
        traces.push_back(Trace("DDDDDMMMIIIIMMMMDDDD", "SSHGNRIVHLQ"));
    }

    AminoAcid abc;
    TraceVector traces;
    size_t length;
};

TEST_F(ProfileBuilderTest, test_build) {
    ProfileBuilder builder(abc);
    Profile profile = builder.build(traces);
    ASSERT_EQ(length, profile.get_length());

    VectorXf col(abc.get_canonical_size());
    col << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
           0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(col.matrix() == profile.get_prob(0).matrix()) << profile.get_prob(0);
    EXPECT_EQ(1, profile.get_eff_num(0));

    col << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0,
           0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0;
    EXPECT_TRUE(col.matrix() == profile.get_prob(7).matrix());
    EXPECT_EQ(2, profile.get_eff_num(7));
}
