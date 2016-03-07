#ifndef _MRFSTATCOMMAND_H_
#define _MRFSTATCOMMAND_H_

#include "command.h"
#include "analyze.h"

class MRFStatCommandLine : public MRFCommandLine {
  public:
    MRFStatCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    struct Option {
        string mrf_filename;
        MRFModelAnalyzer::StatOption stat_opt;
    } opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
    bool parse_mode(char* optarg, StatMode& arg);
    bool parse_corr(char* optarg, StatCorrect& arg);
};

#endif
