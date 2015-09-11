#include <gtest/gtest.h>

#include <sstream>

#include "seqio.h"

class SeqParserTest : public testing::Test {
  public:
    SeqParserTest() {
        SeqRecord r1;
        r1.id   = "d1a6m__";
        r1.desc = "d1a6m__ a.1.1.2 (-) Myoglobin {Sperm whale (Physeter catodon)}";
        r1.seq  = "vlsegewqlvlhvwakveadvaghgqdilirlfkshpetlekfdrfkhlkteaemkased"
                  "lkkhgvtvltalgailkkkghheaelkplaqshatkhkipikylefiseaiihvlhsrhp"
                  "gdfgadaqgamnkalelfrkdiaakykelgy";
        SeqRecord r2;
        r2.id   = "d1mba__";
        r2.desc = "d1mba__ a.1.1.2 (-) Myoglobin {Sea hare (Aplysia limacina)}";
        r2.seq  = "slsaaeadlagkswapvfanknangldflvalfekfpdsanffadfkgksvadikaspkl"
                  "rdvssriftrlnefvnnaanagkmsamlsqfakehvgfgvgsaqfenvrsmfpgfvasva"
                  "appagadaawtklfgliidalkaaga";

        record.push_back(r1);
        record.push_back(r2);
    }

  protected:
    std::vector<SeqRecord> record;
};

TEST_F(SeqParserTest, test_fasta_parser) {
    std::string fasta(
        ">d1a6m__ a.1.1.2 (-) Myoglobin {Sperm whale (Physeter catodon)}\n"
        "vlsegewqlvlhvwakveadvaghgqdilirlfkshpetlekfdrfkhlkteaemkased\n"
        "lkkhgvtvltalgailkkkghheaelkplaqshatkhkipikylefiseaiihvlhsrhp\n"
        "gdfgadaqgamnkalelfrkdiaakykelgy\n"
        "   \n"
        ">d1mba__ a.1.1.2 (-) Myoglobin {Sea hare (Aplysia limacina)}  \n"
        "slsaaeadlagkswapvfanknangldflvalfekfpdsanffadfkgksvadikaspkl  \n"
        "rdvssriftrlnefvnnaanagkmsamlsqfakehvgfgvgsaqfenvrsmfpgfvasva  \n"
        "appagadaawtklfgliidalkaaga\n"
        );

    std::istringstream is(fasta);
    FastaParser parser(is);
    for (size_t i = 0; i < record.size(); ++i) {
        ASSERT_TRUE(parser.has_next());
        SeqRecord r = parser.next();
        EXPECT_EQ(record[i].id, r.id);
        EXPECT_EQ(record[i].desc, r.desc);
        EXPECT_EQ(record[i].seq, r.seq);
    }
    ASSERT_FALSE(parser.has_next());
}
