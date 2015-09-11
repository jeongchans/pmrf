#ifndef _MRFMAINCOMMAND_H_
#define _MRFMAINCOMMAND_H_

#include "command.h"

enum SubCommand { NONE, HELP, BUILD };

class MRFMainCommandLine : public MRFCommandLine {
  public:
    MRFMainCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    SubCommand subcmd;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

#endif
