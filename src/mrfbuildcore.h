#ifndef _MRFBUILDCORE_H_
#define _MRFBUILDCORE_H_

#include "core.h"
#include "mrfbuildcommand.h"

using std::string;

class MRFBuildProcessor : public MRFCmdProcessor {
  public:
    MRFBuildProcessor(int argc, char** argv);
    ~MRFBuildProcessor();

    int build(const string& msa_filename, const string& out_filename);

};

#endif
