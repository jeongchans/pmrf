#ifndef _MRFINFERCORE_H_
#define _MRFINFERCORE_H_

#include "core.h"
#include "mrfinfercommand.h"

using std::string;

class MRFInferProcessor : public MRFCmdProcessor {
  public:
    MRFInferProcessor(int argc, char** argv);
    ~MRFInferProcessor();

    int infer(const string& mrf_filename, const string& seq_filename);

};

#endif
