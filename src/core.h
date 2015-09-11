#ifndef _CORE_H_
#define _CORE_H_

#include <fstream>
#include <string>

#include "protbinfo/seq.h"

#include "command.h"

// AMINO_ACID considers only 20 canonical amino acids
const AminoAcid AMINO_ACID;

// AA_GAP considers 20 canonical amino acids and gap letter
class AminoAcidGap : public Alphabet {
  public:
    AminoAcidGap();

    bool is_canonical_aa(const char& x) const;
};

const AminoAcidGap AA_GAP;

void assert_file_handle(std::ifstream& ifs, const std::string& filename);
void assert_file_handle(std::ofstream& ofs, const std::string& filename);

class MRFCmdProcessor {
  public:
    int run();

  protected:
    MRFCommandLine *cmd_line;
};

#endif
