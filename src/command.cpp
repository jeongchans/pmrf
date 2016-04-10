#include "command.h"

#include <getopt.h>

#include "core.h"

using std::cout;
using std::endl;

bool MRFCommandLine::parse_bool(char* optarg, bool& arg) {
    int n = atoi(optarg);
    if (n == 0) arg = false;
    else if (n == 1) arg = true;
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

bool MRFCommandLine::parse_float(char* optarg, double& arg) {
    arg = atof(optarg);
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
    double regnode_lambda = 0.01;
    double regedge_lambda = 0.2;
    bool regedge_scale = true;
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {"regul", required_argument, 0, 0},
        {"regnode-lambda", required_argument, 0, 0},
        {"regedge-lambda", required_argument, 0, 0},
        {"regedge-scale", required_argument, 0, 0},
        {"msa", required_argument, 0, 0},
        {"edge", required_argument, 0, 0},
        {"delta", required_argument, 0, 0},
        {"seqwt", required_argument, 0, 0},
        {"effnum", required_argument, 0, 0},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "ho:", opts, &opt_idx);
        if (c == -1) break;
        switch (c) {
        case 0:
            switch (opt_idx) {
            case 0:
                show_help();
                exit(0);
            case 1:
                if (parse_regul(optarg, opt.parameterizer_opt.regul)) break;
                else return false;
            case 2:
                if (parse_double(optarg, regnode_lambda)) break;
                else return false;
            case 3:
                if (parse_double(optarg, regedge_lambda)) break;
                else return false;
            case 4:
                if (parse_bool(optarg, regedge_scale)) break;
                else return false;
            case 5:
                if (parse_msa_fmt(optarg, opt.msa_fmt)) break;
                else return false;
            case 6:
                if (parse_str(optarg, opt.eidx_filename)) break;
                else return false;
            case 7:
                if (parse_float(optarg, opt.optim_opt.delta)) break;
                else return false;
            case 8:
                if (parse_int(optarg, opt.msa_analyzer_opt.seq_wt)) break;
                else return false;
            case 9:
                if (parse_int(optarg, opt.msa_analyzer_opt.eff_num)) break;
                else return false;
            }
            break;
        case 'h':
            show_help();
            exit(0);
        case 'o':
            opt.out_filename = optarg;
            break;
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
    if (opt.parameterizer_opt.regul == RegulMethod::RegulMethod::L2) {
        opt.parameterizer_opt.l2_opt.lambda1 = regnode_lambda;
        opt.parameterizer_opt.l2_opt.lambda2 = regedge_lambda;
        opt.parameterizer_opt.l2_opt.sc = regedge_scale;
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
    int d = atoi(optarg);
    if (d == 0) arg = RegulMethod::RegulMethod::NONE;
    else if (d == 1) arg = RegulMethod::RegulMethod::L2;
    else {
        error_message = string("Unknown node regularization option: ") + *optarg;
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::parse_double(char* optarg, double& arg) {
    arg = atof(optarg);
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
