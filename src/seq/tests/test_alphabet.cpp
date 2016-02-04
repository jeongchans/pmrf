#include <gtest/gtest.h>
#include <blitz/array.h>

#include <string>
#include <vector>

#include "alphabet.h"

class AlphabetTest : public testing::Test {
  protected:
    virtual void SetUp() {
        canonical_sym = "ABC";
        gap_sym = "-";
        degenerate_sym = "PR";
        unknown_sym = "X";
        none_sym = "*";
        missing_sym = "~";
        undefined_sym = "LMN@#$%^&()<>?";
        abc = Alphabet(canonical_sym.c_str(),
                       gap_sym.c_str(),
                       degenerate_sym.c_str(),
                       unknown_sym.c_str(),
                       none_sym.c_str(),
                       missing_sym.c_str());
        abc.set_degeneracy('P', 'A');
        abc.set_degeneracy('P', 'B');
        abc.set_degeneracy('R', 'C');
    }

    Alphabet abc;
    std::string canonical_sym;
    std::string gap_sym;
    std::string degenerate_sym;
    std::string unknown_sym;
    std::string none_sym;
    std::string missing_sym;
    std::string undefined_sym;
};

TEST_F(AlphabetTest, test_get_canonical) {
    EXPECT_EQ(canonical_sym, abc.get_canonical());
}

TEST_F(AlphabetTest, test_get_gap) {
    EXPECT_EQ(gap_sym, abc.get_gap());
}

TEST_F(AlphabetTest, test_get_unknown) {
    EXPECT_EQ(unknown_sym, abc.get_unknown());
}

TEST_F(AlphabetTest, test_is_canonical) {
    std::string::iterator pos;
    for (pos = canonical_sym.begin(); pos != canonical_sym.end(); ++pos) {
        EXPECT_TRUE(abc.is_canonical(*pos));
    }
    std::string s = gap_sym + degenerate_sym + unknown_sym + none_sym + missing_sym + undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_canonical(*pos));
    }
}

TEST_F(AlphabetTest, test_is_gap) {
    std::string::iterator pos;
    for (pos = gap_sym.begin(); pos != gap_sym.end(); ++pos) {
        EXPECT_TRUE(abc.is_gap(*pos));
    }
    std::string s = canonical_sym + degenerate_sym + unknown_sym + none_sym + missing_sym + undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_gap(*pos));
    }
}

TEST_F(AlphabetTest, test_is_degenerate) {
    std::string::iterator pos;
    for (pos = degenerate_sym.begin(); pos != degenerate_sym.end(); ++pos) {
        EXPECT_TRUE(abc.is_degenerate(*pos));
    }
    std::string s = canonical_sym + gap_sym + unknown_sym + none_sym + missing_sym + undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_degenerate(*pos));
    }
}

TEST_F(AlphabetTest, test_is_unknown) {
    std::string::iterator pos;
    for (pos = unknown_sym.begin(); pos != unknown_sym.end(); ++pos) {
        EXPECT_TRUE(abc.is_unknown(*pos));
    }
    std::string s = canonical_sym + gap_sym + degenerate_sym + none_sym + missing_sym + undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_unknown(*pos));
    }
}

TEST_F(AlphabetTest, test_is_none) {
    std::string::iterator pos;
    for (pos = none_sym.begin(); pos != none_sym.end(); ++pos) {
        EXPECT_TRUE(abc.is_none(*pos));
    }
    std::string s = canonical_sym + gap_sym + degenerate_sym + unknown_sym + missing_sym + undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_none(*pos));
    }
}

TEST_F(AlphabetTest, test_is_missing) {
    std::string::iterator pos;
    for (pos = missing_sym.begin(); pos != missing_sym.end(); ++pos) {
        EXPECT_TRUE(abc.is_missing(*pos));
    }
    std::string s = canonical_sym + gap_sym + degenerate_sym + unknown_sym + none_sym + undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_missing(*pos));
    }
}

TEST_F(AlphabetTest, test_get_idx) {
    std::string s = canonical_sym + gap_sym + degenerate_sym + unknown_sym + none_sym + missing_sym;
    for (int i = 0; i < (int) s.size(); ++i) {
        ASSERT_EQ(i, abc.get_idx(s[i]));
    }
}

TEST_F(AlphabetTest, test_case_insensitive) {
    Alphabet new_abc(canonical_sym.c_str(),
                     gap_sym.c_str(),
                     degenerate_sym.c_str(),
                     unknown_sym.c_str(),
                     none_sym.c_str(),
                     missing_sym.c_str(),
                     true);
    std::string s = canonical_sym + gap_sym + degenerate_sym + unknown_sym + none_sym + missing_sym;
    std::string::iterator pos;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_EQ(new_abc.is_canonical(toupper(*pos)), new_abc.is_canonical(tolower(*pos)));
        EXPECT_EQ(new_abc.is_gap(toupper(*pos)), new_abc.is_gap(tolower(*pos)));
        EXPECT_EQ(new_abc.is_degenerate(toupper(*pos)), new_abc.is_degenerate(tolower(*pos)));
        EXPECT_EQ(new_abc.is_unknown(toupper(*pos)), new_abc.is_unknown(tolower(*pos)));
        EXPECT_EQ(new_abc.is_none(toupper(*pos)), new_abc.is_none(tolower(*pos)));
        EXPECT_EQ(new_abc.is_missing(toupper(*pos)), new_abc.is_missing(tolower(*pos)));
        ASSERT_EQ(new_abc.get_idx(toupper(*pos)), new_abc.get_idx(tolower(*pos)));
    }
}

TEST_F(AlphabetTest, test_is_valid) {
    std::string::iterator pos;
    std::string s = canonical_sym + gap_sym + degenerate_sym + unknown_sym + none_sym + missing_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(abc.is_valid(*pos));
    }
    s = undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(abc.is_valid(*pos));
    }
}

TEST_F(AlphabetTest, test_get_canonical_size) {
    EXPECT_EQ(canonical_sym.size(), abc.get_canonical_size());
}

TEST_F(AlphabetTest, test_get_degeneracy) {
    EXPECT_EQ(std::string("AB"), abc.get_degeneracy('P'));
    EXPECT_EQ(std::string("C"), abc.get_degeneracy('R'));

    FloatType w;
    EXPECT_EQ(std::string("AB"), abc.get_degeneracy('P', &w));
    EXPECT_FLOAT_EQ(0.5, w);
    EXPECT_EQ(std::string("C"), abc.get_degeneracy('R', &w));
    EXPECT_FLOAT_EQ(1.0, w);
    EXPECT_EQ(std::string(""), abc.get_degeneracy('X', &w));
    EXPECT_FLOAT_EQ(0.0, w);
}

TEST_F(AlphabetTest, test_get_count) {
    ASSERT_EQ(canonical_sym.size(), (size_t)abc.get_count('A').size());
    blitz::Array<double, 1> cnt(3);
    cnt = 1, 0, 0;
    EXPECT_TRUE(all(cnt == abc.get_count('A')));
    cnt = 0.5, 0.5, 0;
    EXPECT_TRUE(all(cnt == abc.get_count('P')));
    cnt = 0, 0, 1;
    EXPECT_TRUE(all(cnt == abc.get_count('R')));
    cnt = 0, 0, 0;
    EXPECT_TRUE(all(cnt == abc.get_count('X')));
    EXPECT_TRUE(all(cnt == abc.get_count('-')));
    EXPECT_TRUE(all(cnt == abc.get_count('*')));
    EXPECT_TRUE(all(cnt == abc.get_count('~')));
}

class AminoAcidTest : public testing::Test {
  protected:
    virtual void SetUp() {
        valid_sym = "ACDEFGHIKLMNPQRSTVWY-BJZOUX*~";
        undefined_sym = "0123456789!@#$%^&()<>?";
        amino = AminoAcid();
    }

    AminoAcid amino;
    std::string valid_sym;
    std::string undefined_sym;

};

TEST_F(AminoAcidTest, test_is_xxx) {
    std::string::iterator pos;
    std::string s;
    s = valid_sym.substr(0, 20);
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_canonical(*pos));
    }
    s = valid_sym.substr(20, 1);
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_gap(*pos));
    }
    s = valid_sym.substr(21, 5);
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_degenerate(*pos));
    }
    s = valid_sym.substr(26, 1);
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_unknown(*pos));
    }
    s = valid_sym.substr(27, 1);
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_none(*pos));
    }
    s = valid_sym.substr(28, 1);
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_missing(*pos));
    }
    s = valid_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_TRUE(amino.is_valid(*pos));
    }
    s = undefined_sym;
    for (pos = s.begin(); pos != s.end(); ++pos) {
        EXPECT_FALSE(amino.is_valid(*pos));
    }
}

TEST_F(AminoAcidTest, test_get_idx) {
    for (int i = 0; i < (int) valid_sym.size(); ++i) {
        ASSERT_EQ(i, amino.get_idx(valid_sym[i]));
    }
}

TEST_F(AminoAcidTest, test_get_count) {
    ASSERT_EQ(20, amino.get_count('C').size());
    EXPECT_EQ(1, sum(amino.get_count('C')));
    EXPECT_EQ(1, amino.get_count('C')(1));
    EXPECT_EQ(1, sum(amino.get_count('B')));
    EXPECT_EQ(0.5, amino.get_count('B')(2));
    EXPECT_EQ(0.5, amino.get_count('B')(11));
}
