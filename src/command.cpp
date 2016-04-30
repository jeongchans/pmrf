#include "command.h"

#include <getopt.h>

#include "core.h"

using std::cout;
using std::endl;

bool MRFCommandLine::parse_bool(char* optarg, bool& arg) {
    string val = string(optarg);
    if (val == "yes") arg = true;
    else if (val == "no") arg = false;
    else {
        error_message = "Not acceptable: " + string(optarg);
        return false;
    }
    return true;
}

bool MRFCommandLine::parse_int(char* optarg, int& arg) {
    arg = atoi(optarg);
    return true;
}

bool MRFCommandLine::parse_float(char* optarg, float& arg) {
    arg = atof(optarg);
    return true;
}

bool MRFCommandLine::parse_str(char* optarg, string& arg) {
    arg = string(optarg);
    return true;
}

bool MRFCommandLine::set_opt_err_msg(const string& opt, const char* optarg) {
    error_message = string("Not acceptable option: '") + opt + " " + optarg + "'\n";
    return false;
}

/**
   @class MRFMainCommandLine
 */

MRFMainCommandLine::MRFMainCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    subcmd = NONE;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFMainCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFMainProcessor*)processor)->run_mrf_cmd(subcmd);
}

void MRFMainCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " <command> [<args>]" << endl
         << endl
         << "The following commands are available for MRF modeling and varisous applications:" << endl
         << endl
         << Main::command_list_message
         << endl;
}

void MRFMainCommandLine::show_version() {
    cout << PROGNAME << " version " << VERSION << endl;
}

bool MRFMainCommandLine::parse_command_line(int argc, char** argv) {
    if (argc < 2) {
        show_help();
        exit(0);
    }
    if (argv[1][0] != '-') {
        string cmd = argv[1];
        if (cmd == "help") subcmd = HELP;
        else if (cmd == "build") subcmd = BUILD;
        else if (cmd == "infer") subcmd = INFER;
        else if (cmd == "stat") subcmd = STAT;
        else if (cmd == "show") subcmd = SHOW;
        else {
            error_message = string("Unknown command: ") + cmd + "\n";
            return false;
        }
        return true;
    }
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {"version", 0, 0, 0},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        switch (c) {
        case 0:
            switch (opt_idx) {
            case 0:
                show_help();
                exit(0);
            case 1:
                show_version();
                exit(0);
            }
            break;
        case 'h':
            show_help();
            exit(0);
        default:
            return false;
        }
    }
    return true;
}

/**
   @class MRFBuildCommandLine
 */

MRFBuildCommandLine::MRFBuildCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    opt.out_filename = "";
    opt.msa_fmt = AFASTA;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFBuildCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFBuildProcessor*)processor)->build();
}

void MRFBuildCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " build <msa_file> [options]" << endl
         << endl
         << Build::option_message
         << endl;
}

bool MRFBuildCommandLine::parse_command_line(int argc, char** argv) {
    string sval;
    float fval;
    int ival;
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        /* official options */
        {"help",            no_argument,        0, 'h'},
        {"msa",             required_argument,  0, 100},
        {"edge",            required_argument,  0, 'e'},
        {"output",          required_argument,  0, 'o'},
        {"seqwt",           required_argument,  0, 300},
//        {"neff",            required_argument,  0, 0},
        {"clstr-maxidt",    required_argument,  0, 310},
        {"symmetric",       no_argument,        0, 400},
        {"regul",           required_argument,  0, 401},
        {"regul-vl",        required_argument,  0, 410},
        {"regul-wl",        required_argument,  0, 411},
        {"lbfgs-corr",      required_argument,  0, 510},
        {"lbfgs-epsilon",   required_argument,  0, 511},
        {"lbfgs-delta",     required_argument,  0, 512},
        {"lbfgs-maxiter",   required_argument,  0, 513},
        /* experimental options */
//        {"termi-maxgap",    required_argument,  0, -100},
        {"regw-lambda-max", required_argument,  0, -200},
        {"regw-lambda-min", required_argument,  0, -201},
        {"regw-lambda-sc",  required_argument,  0, -202},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "e:o:h", opts, &opt_idx);
        if (c == -1) break;
        /* official options */
        else if (c == 100) validity = parse_str(optarg, sval) && set_msa_fmt(opt.msa_fmt, sval);
        else if (c == 'e') validity = parse_str(optarg, opt.eidx_filename);
        else if (c == 'o') validity = parse_str(optarg, opt.out_filename);
        else if (c == 300) validity = parse_str(optarg, sval) && set_seq_wt(opt.msa_analyzer_opt.seq_wt, sval);
        else if (c == 310) validity = parse_float(optarg, fval) && fval >= 0. && fval <= 1. && set_val<float>(opt.msa_analyzer_opt.clstr_maxidt, fval);
        else if (c == 400) validity = set_val<bool>(opt.parameterizer_opt.asymmetric, !opt.parameterizer_opt.asymmetric);
        else if (c == 401) validity = parse_str(optarg, sval) && set_regul(opt.parameterizer_opt.regul, sval);
        else if (c == 410) validity = parse_float(optarg, fval) && fval >= 0. && set_val<float>(opt.parameterizer_opt.regnode_lambda, fval);
        else if (c == 411) validity = parse_float(optarg, fval) && fval >= 0. && set_val<float>(opt.parameterizer_opt.regedge_lambda, fval);
        else if (c == 510) validity = parse_int(optarg, ival) && ival >= 3 && set_val<int>(opt.optim_opt.corr, ival);
        else if (c == 511) validity = parse_float(optarg, fval) && fval > 0. && set_val<float>(opt.optim_opt.epsilon, fval);
        else if (c == 512) validity = parse_float(optarg, fval) && fval > 0. && set_val<float>(opt.optim_opt.delta, fval);
        else if (c == 513) validity = parse_int(optarg, ival) && ival > 0 && set_val<int>(opt.optim_opt.max_iter, ival);
        else if (c == 'h') { show_help(); exit(0); }
        /* experimental options */
//        else if (c == -100) validity = parse_float(optarg, opt.msa_analyzer_opt.termi_maxgapperc);
        else if (c == -200) validity = parse_float(optarg, opt.parameterizer_opt.regedge_lambda_max);
        else if (c == -201) validity = parse_float(optarg, opt.parameterizer_opt.regedge_lambda_min);
        else if (c == -202) validity = parse_float(optarg, opt.parameterizer_opt.regedge_lambda_sc);
        else validity = false;
        if (!validity) return set_opt_err_msg("--" + string(opts[opt_idx].name), optarg);
    }
    if (optind != argc - 1) { show_help(); exit(EXIT_FAILURE); }
    opt.msa_filename = argv[optind++];
    return true;
}

bool MRFBuildCommandLine::set_msa_fmt(MSAFormat& arg, const string& val) {
    if (val == "fasta") arg = AFASTA;
    else if (val == "a3m") arg = A3M;
    else return false;
    return true;
}

bool MRFBuildCommandLine::set_seq_wt(MSAProcOption::SeqWeight& arg, const string& val) {
    if (val == "no") arg = MSAProcOption::SW_NO;
    else if (val == "clstr") arg = MSAProcOption::SW_CLSTR;
    else if (val == "pb") arg = MSAProcOption::SW_PB;
    else return false;
    return true;
}

//bool MRFBuildCommandLine::parse_eff_num(char* optarg, MSAProcOption::EffSeqNum& arg) {
//    string val = string(optarg);
//    if (val == "no") arg = MSAProcOption::NEFF_NO;
//    else if (val == "clstr") arg = MSAProcOption::NEFF_CLSTR;
//    else {
//        error_message = "Unsupported option for effective number estimation: " + val;
//        return false;
//    }
//    return true;
//}

bool MRFBuildCommandLine::set_regul(RegulMethod::RegulMethod& arg, const string& val) {
    if (val == "no") arg = RegulMethod::REGUL_NONE;
    else if (val == "l2") arg = RegulMethod::REGUL_L2;
    else return false;
    return true;
}

/**
   @class MRFStatCommandLine
 */

MRFStatCommandLine::MRFStatCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFStatCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFStatProcessor*)processor)->stat();
}

void MRFStatCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " stat <mrf_file> [options]" << endl
         << endl
         << Stat::option_message
         << endl;
}

bool MRFStatCommandLine::parse_command_line(int argc, char** argv) {
    string sval;
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help",            no_argument,        0, 'h'},
        {"mode",            required_argument,  0, 100},
        {"corr",            required_argument,  0, 300},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        else if (c == 100) validity = parse_str(optarg, sval) && set_mode(opt.mode, sval);
        else if (c == 300) validity = parse_str(optarg, sval) && set_corr(opt.corr, sval);
        else if (c == 'h') { show_help(); exit(0); }
        else validity = false;
        if (!validity) return set_opt_err_msg("--" + string(opts[opt_idx].name), optarg);
    }
    if (optind != argc - 1) { show_help(); exit(EXIT_FAILURE); }
    opt.mrf_filename = argv[optind++];
    return true;
}

bool MRFStatCommandLine::set_mode(Stat::Mode& arg, const string& val) {
    if (val == "pair") arg = Stat::MODE_PAIR;
    else if (val == "pos") arg = Stat::MODE_POS;
    else return false;
    return true;
}

bool MRFStatCommandLine::set_corr(Stat::Correct& arg, const string& val) {
    if (val == "no") arg = Stat::CORR_NONE;
    else if (val == "apc") arg = Stat::CORR_APC;
    else if (val == "ncps") arg = Stat::CORR_NCPS;
    else return false;
    return true;
}

/**
   @class MRFInferCommandLine
 */

MRFInferCommandLine::MRFInferCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFInferCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFInferProcessor*)processor)->infer();
}

void MRFInferCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " infer <mrf_file> <seq_file> [options]" << endl
         << endl
         << Infer::option_message
         << endl;
}

bool MRFInferCommandLine::parse_command_line(int argc, char** argv) {
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        switch (c) {
        case 0:
            switch (opt_idx) {
            case 0:
                show_help();
                exit(0);
            }
            break;
        case 'h':
            show_help();
            exit(0);
        default:
            return false;
        }
    }
    if (optind < argc) {
        opt.mrf_filename = argv[optind++];
        opt.seq_filename = argv[optind];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    return true;
}

/**
   @class MRFShowCommandLine
 */

MRFShowCommandLine::MRFShowCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFShowCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFShowProcessor*)processor)->show();
}

void MRFShowCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " show <mrf_file> [options]" << endl
         << endl
         << Show::option_message
         << endl;
}

bool MRFShowCommandLine::parse_command_line(int argc, char** argv) {
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        switch (c) {
        case 0:
            switch (opt_idx) {
            case 0:
                show_help();
                exit(0);
            }
            break;
        case 'h':
            show_help();
            exit(0);
        default:
            return false;
        }
    }
    if (optind < argc) {
        opt.mrf_filename = argv[optind++];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    return true;
}
