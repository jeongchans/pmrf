#ifndef _MRFINFERCOMMAND_H_
#define _MRFINFERCOMMAND_H_

#include "command.h"
#include "analyze.h"

class MRFInferCommandLine : public MRFCommandLine {
  public:
    MRFInferCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    struct Option {
        string mrf_filename;
        string seq_filename;
    } opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

#endif

