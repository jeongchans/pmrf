#ifndef _COMMAND_H_
#define _COMMAND_H_

#include <string>
#include <iostream>
#include <cstring>

#include "option.h"

#define QUOTE(name)     #name
#define PROGNAME        QUOTE(pmrf)
#define VERSION         QUOTE(0.2.0)

using std::string;

class MRFCmdProcessor;

class MRFCommandLine {
  public:
    MRFCommandLine(int argc, char** argv) : validity(false), error_message(""), argc(argc) {
        this->argv = (char**)malloc(sizeof(char*) * argc);
        for (int i = 0; i < argc; ++i) {
            this->argv[i] = (char*)malloc(strlen(argv[i]) + 1);
            strcpy(this->argv[i], argv[i]);
        }
    }

    virtual ~MRFCommandLine() {
        for (int i = 0; i < argc; ++i) free(argv[i]);
        free(argv);
    }

    virtual int process_command(MRFCmdProcessor *processor) = 0;
    virtual void show_help() = 0;

    bool is_valid() const { return validity; }
    void show_error() { std::cerr << error_message << std::endl; }

  protected:
    bool validity;
    string error_message;
    int argc;
    char **argv;

    virtual bool parse_command_line(int argc, char** argv) = 0;

    bool parse_bool(char* optarg, bool& arg);
    bool parse_int(char* optarg, int& arg);
    bool parse_float(char* optarg, double& arg);
    bool parse_str(char* optarg, string& arg);
};

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

class MRFBuildCommandLine : public MRFCommandLine {
  public:
    MRFBuildCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    Build::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
    bool parse_msa_fmt(char* optarg, MSAFormat& arg);
    bool parse_regul(char* optarg, RegulMethod::RegulMethod& arg);
    bool parse_double(char* optarg, double& arg);
};

class MRFStatCommandLine : public MRFCommandLine {
  public:
    MRFStatCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    Stat::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
    bool parse_mode(char* optarg, Stat::Mode& arg);
    bool parse_corr(char* optarg, Stat::Correct& arg);
};

class MRFInferCommandLine : public MRFCommandLine {
  public:
    MRFInferCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);
    virtual void show_help();

    Infer::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

#endif
