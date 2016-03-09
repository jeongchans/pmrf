#ifndef _MRFMAINCOMMAND_H_
#define _MRFMAINCOMMAND_H_

#include "command.h"

enum SubCommand { NONE, HELP, BUILD, INFER, STAT };

class MRFMainCommandLine : public MRFCommandLine {
  public:
    MRFMainCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    void show_version();

    SubCommand subcmd;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

#endif
