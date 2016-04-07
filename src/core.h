#ifndef _CORE_H_
#define _CORE_H_

#include <fstream>
#include <string>

#include "seq/alphabet.h"

#include "command.h"

using std::string;

const AminoAcid AA("-", false, true);
//const AminoAcid AA("=-^", false, true);
//const char GAP_OPEN_SYMBOL = '=';
//const char GAP_EXT_SYMBOL = '-';
//const char GAP_UNALI_SYMBOL = '^';

void assert_file_handle(std::ifstream& ifs, const std::string& filename);
void assert_file_handle(std::ofstream& ofs, const std::string& filename);

class MRFCommandLine;

class MRFCmdProcessor {
  public:
    int run();

  protected:
    MRFCommandLine *cmd_line;
};

class MRFMainProcessor : public MRFCmdProcessor {
  public:
    MRFMainProcessor(int argc, char** argv);
    ~MRFMainProcessor();

    int run_mrf_cmd(const SubCommand& cmd);

  private:
    int argc;
    char** argv;
};

class MRFBuildProcessor : public MRFCmdProcessor {
  public:
    MRFBuildProcessor(int argc, char** argv);
    ~MRFBuildProcessor();

    int build(const string& msa_filename, const string& out_filename);

};

class MRFStatProcessor : public MRFCmdProcessor {
  public:
    MRFStatProcessor(int argc, char** argv);
    ~MRFStatProcessor();

    int stat(const string& mrf_filename);

};

class MRFInferProcessor : public MRFCmdProcessor {
  public:
    MRFInferProcessor(int argc, char** argv);
    ~MRFInferProcessor();

    int infer(const string& mrf_filename, const string& seq_filename);

};

#endif
