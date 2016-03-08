#ifndef _MRFSTATCORE_H_
#define _MRFSTATCORE_H_

#include "core.h"
#include "mrfstatcommand.h"

using std::string;

class MRFStatProcessor : public MRFCmdProcessor {
  public:
    MRFStatProcessor(int argc, char** argv);
    ~MRFStatProcessor();

    int stat(const string& mrf_filename);

};

#endif
