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

    MRFCommandLine(int argc, char** argv, const string& usage_message, const string& opt_message);
    virtual ~MRFCommandLine();

    virtual int process_command(MRFCmdProcessor *processor) = 0;

    void show_usage() { std::cout << "Usage: " << PROGNAME << " " << usage_message << std::endl; }
    void show_help() { show_usage(); std::cout << std::endl << opt_message << std::endl; }
    bool is_valid() const { return validity; }
    void show_error() { std::cerr << error_message << std::endl; }

  protected:
    int argc;
    char **argv;
    const string usage_message;
    const string opt_message;
    bool validity;
    string error_message;

    virtual bool parse_command_line(int argc, char** argv) = 0;

    bool parse_int(char* optarg, int& arg);
    bool parse_float(char* optarg, float& arg);
    bool parse_str(char* optarg, string& arg);
    template <typename T> bool set_val(T& arg, const T& val) { arg = val; return true; }

    bool set_opt_err_msg(const string& opt, const char* optarg);
};

class MRFMainCommandLine : public MRFCommandLine {
  public:
    MRFMainCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);

    void show_version() { std::cout << PROGNAME << " version " << VERSION << std::endl; }

    SubCommand subcmd;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

class MRFBuildCommandLine : public MRFCommandLine {
  public:
    MRFBuildCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);

    Build::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
    bool set_msa_fmt(MSAFormat& arg, const string& val);
    bool set_seq_wt(MSAProcOption::SeqWeight& arg, const string& val);
    bool set_eff_num(MSAProcOption::EffSeqNum& arg, const string& val);
    bool set_regul(RegulMethod::RegulMethod& arg, const string& val);
};

class MRFStatCommandLine : public MRFCommandLine {
  public:
    MRFStatCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);

    Stat::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
    bool set_mode(Stat::Mode& arg, const string& val);
    bool set_corr(Stat::Correct& arg, const string& val);
};

class MRFInferCommandLine : public MRFCommandLine {
  public:
    MRFInferCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);

    Infer::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

class MRFShowCommandLine : public MRFCommandLine {
  public:
    MRFShowCommandLine(int argc, char** argv);

    virtual int process_command(MRFCmdProcessor *processor);

    Show::Option opt;

  protected:
    virtual bool parse_command_line(int argc, char** argv);
};

#endif
