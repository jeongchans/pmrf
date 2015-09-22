#ifndef _CORE_H_
#define _CORE_H_

#include <fstream>
#include <string>

#include "seq/alphabet.h"

#include "command.h"

//const AminoAcid AA("-", false, true);
const AminoAcid AA("=-^", false, true);

void assert_file_handle(std::ifstream& ifs, const std::string& filename);
void assert_file_handle(std::ofstream& ofs, const std::string& filename);

class MRFCmdProcessor {
  public:
    int run();

  protected:
    MRFCommandLine *cmd_line;
};

#endif
