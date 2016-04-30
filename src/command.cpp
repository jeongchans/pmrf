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

bool MRFCommandLine::parse_float(char* optarg, float& arg, const float& minval, const float& maxval) {
    arg = atof(optarg);
    if (arg < minval || arg > maxval) {
        error_message = "Not acceptable: " + string(optarg) + "\n" +
                        "The value should be in the range from " + std::to_string(minval) + " to " + std::to_string(maxval);
        return false;
    }
    return true;
}

bool MRFCommandLine::parse_str(char* optarg, string& arg) {
    arg = string(optarg);
    return true;
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
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        /* official options */
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
//        {"regw-sc-deg",     required_argument,  0, 0},
//        {"regw-sc-neff",    required_argument,  0, 0},
        {"lbfgs-corr",      required_argument,  0, 510},
        {"lbfgs-epsilon",   required_argument,  0, 511},
        {"lbfgs-delta",     required_argument,  0, 512},
        {"lbfgs-maxiter",   required_argument,  0, 513},
        {"help",            no_argument,        0, 'h'},
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
        switch (c) {
        /* official options */
        case 100:
            if (parse_msa_fmt(optarg, opt.msa_fmt)) break;
            else return false;
        case 'e':
            if (parse_str(optarg, opt.eidx_filename)) break;
            else return false;
        case 'o':
            opt.out_filename = optarg;
            break;
        case 300:
            if (parse_seq_wt(optarg, opt.msa_analyzer_opt.seq_wt)) break;
            else return false;
        case 310:
            if (parse_float(optarg, opt.msa_analyzer_opt.clstr_maxidt, 0., 1.)) break;
            else return false;
//        case 9:
//            if (parse_eff_num(optarg, opt.msa_analyzer_opt.eff_num)) break;
//            else return false;
        case 400:
            opt.parameterizer_opt.asymmetric = !opt.parameterizer_opt.asymmetric;
            break;
        case 401:
            if (parse_regul(optarg, opt.parameterizer_opt.regul)) break;
            else return false;
        case 410:
            if (parse_float(optarg, opt.parameterizer_opt.regnode_lambda)) break;
            else return false;
        case 411:
            if (parse_float(optarg, opt.parameterizer_opt.regedge_lambda)) break;
            else return false;
//        case 10:
//            //if (parse_bool(optarg, opt.parameterizer_opt.regedge_sc_neff)) break;
//            //else return false;
        case 510:
            if (parse_int(optarg, opt.optim_opt.corr)) break;
            else return false;
        case 511:
            if (parse_float(optarg, opt.optim_opt.epsilon)) break;
            else return false;
        case 512:
            if (parse_float(optarg, opt.optim_opt.delta)) break;
            else return false;
        case 513:
            if (parse_int(optarg, opt.optim_opt.max_iter)) break;
            else return false;
        case 'h':
            show_help();
            exit(0);
        /* experimental options */
//        case -100: parse_float(optarg, opt.msa_analyzer_opt.termi_maxgapperc); break;
        case -200: parse_float(optarg, opt.parameterizer_opt.regedge_lambda_max); break;
        case -201: parse_float(optarg, opt.parameterizer_opt.regedge_lambda_min); break;
        case -202: parse_float(optarg, opt.parameterizer_opt.regedge_lambda_sc); break;
        default:
            return false;
        }
    }
    if (optind < argc) {
        opt.msa_filename = argv[optind];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::parse_msa_fmt(char* optarg, MSAFormat& arg) {
    string fmt = string(optarg);
    if (fmt == "fasta") arg = AFASTA;
    else if (fmt == "a3m") arg = A3M;
    else {
        error_message = "Unsupported MSA format: " + fmt;
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::parse_regul(char* optarg, RegulMethod::RegulMethod& arg) {
    string val = string(optarg);
    if (val == "no") arg = RegulMethod::REGUL_NONE;
    else if (val == "l2") arg = RegulMethod::REGUL_L2;
    else {
        error_message = "Unknown option for regularization: " + val;
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::parse_seq_wt(char* optarg, MSAProcOption::SeqWeight& arg) {
    string val = string(optarg);
    if (val == "no") arg = MSAProcOption::SW_NO;
    else if (val == "clstr") arg = MSAProcOption::SW_CLSTR;
    else if (val == "pb") arg = MSAProcOption::SW_PB;
    else {
        error_message = "Unsupported option for sequence weighting: " + val;
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::parse_eff_num(char* optarg, MSAProcOption::EffSeqNum& arg) {
    string val = string(optarg);
    if (val == "no") arg = MSAProcOption::NEFF_NO;
    else if (val == "clstr") arg = MSAProcOption::NEFF_CLSTR;
    else {
        error_message = "Unsupported option for effective number estimation: " + val;
        return false;
    }
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
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {"mode", required_argument, 0, 0},
        {"corr", required_argument, 0, 0},
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
                if (parse_mode(optarg, opt.mode)) break;
                return false;
            case 2:
                if (parse_corr(optarg, opt.corr)) break;
                return false;
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
        opt.mrf_filename = argv[optind];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    return true;
}

bool MRFStatCommandLine::parse_mode(char* optarg, Stat::Mode& arg) {
    string mode = string(optarg);
    if (mode == "pair") arg = Stat::MODE_PAIR;
    else if (mode == "pos") arg = Stat::MODE_POS;
    else {
        error_message = "Unsupported mode: " + mode;
        return false;
    }
    return true;
}

bool MRFStatCommandLine::parse_corr(char* optarg, Stat::Correct& arg) {
    int d;
    if (!parse_int(optarg, d)) return false;
    if (d == 0) arg = Stat::CORR_NONE;
    else if (d == 1) arg = Stat::CORR_APC;
    else if (d == 2) arg = Stat::CORR_NCPS;
    else {
        error_message = "Unsupported correction: " + string(optarg);
        return false;
    }
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
