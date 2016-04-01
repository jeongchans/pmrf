#include <gtest/gtest.h>

#include "bgfreq.h"

class BgFreqTest : public testing::Test {
  public:
    AminoAcid abc;
};

TEST_F(BgFreqTest, test_bgfreq) {
    std::string s = abc.get_canonical();
    std::map<char, double> freq;
    VectorXf f(s.size());
    for (int i = 0; i < (int)s.size(); ++i) {
        double val = i + 1;
        freq.insert(make_pair(s[i], val));
        f(i) = val;
    }
    for (std::map<char, double>::iterator pos = freq.begin(); pos != freq.end(); ++pos) pos->second /= f.sum();
    f /= f.sum();

    BgFreq bgfreq(freq);
    EXPECT_TRUE(f.matrix() == bgfreq.get_array(abc).matrix());

    for (int i = 0; i < (int)s.size(); ++i) {
        EXPECT_FLOAT_EQ(f(i), bgfreq.get_value(s[i]));
    }
}

TEST_F(BgFreqTest, test_robinson) {
    RobinsonBgFreq robinson;
    VectorXf v(20);
    v << 0.07805,  // A
         0.01925,  // C
         0.05364,  // D
         0.06295,  // E
         0.03856,  // F
         0.07377,  // G
         0.02199,  // H
         0.05142,  // I
         0.05744,  // K
         0.09019,  // L
         0.02243,  // M
         0.04487,  // N
         0.05203,  // P
         0.04264,  // Q
         0.05129,  // R
         0.07120,  // S
         0.05841,  // T
         0.06441,  // V
         0.01330,  // W
         0.03216;  // Y
    EXPECT_TRUE(v.matrix() == robinson.get_array(abc).matrix());
    EXPECT_EQ(0.02199, robinson.get_value('H'));
    EXPECT_EQ(0.07120, robinson.get_value('S'));
}

TEST_F(BgFreqTest, test_altschul) {
    AltschulBgFreq altschul;
    VectorXf v(20);
    v << 0.08100,  // A
         0.01500,  // C
         0.05400,  // D
         0.06100,  // E
         0.04000,  // F
         0.06800,  // G
         0.02200,  // H
         0.05700,  // I
         0.05600,  // K
         0.09300,  // L
         0.02500,  // M
         0.04500,  // N
         0.04900,  // P
         0.03900,  // Q
         0.05700,  // R
         0.06800,  // S
         0.05800,  // T
         0.06700,  // V
         0.01300,  // W
         0.03200;  // Y
    EXPECT_TRUE(v.matrix() == altschul.get_array(abc).matrix());
    EXPECT_EQ(0.02200, altschul.get_value('H'));
    EXPECT_EQ(0.06800, altschul.get_value('S'));
}
