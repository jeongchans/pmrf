#ifndef _MRFBUILDCOMMAND_H_
#define _MRFBUILDCOMMAND_H_

#include "command.h"
#include "build.h"

class MRFBuildCommandLine : public MRFCommandLine {
  public:
    MRFBuildCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    struct Option {
        string msa_filename;
        string out_filename;
        MRFModelBuilder::Option build_opt;
    } opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
    bool parse_msa_fmt(char* optarg, MSAFormat& arg);
    bool parse_regul(char* optarg, RegulMethod::RegulMethod& arg);
    bool parse_double(char* optarg, double& arg);
};

#endif
