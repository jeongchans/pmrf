#ifndef _MRFMAINCORE_H_
#define _MRFMAINCORE_H_

#include "core.h"
#include "mrfmaincommand.h"

using std::string;

class MRFMainProcessor : public MRFCmdProcessor {
  public:
    MRFMainProcessor(int argc, char** argv);
    ~MRFMainProcessor();

    int run_mrf_cmd(const SubCommand& cmd);

  private:
    int argc;
    char** argv;
};

#endif
