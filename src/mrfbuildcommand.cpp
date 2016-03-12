#include "mrfbuildcommand.h"
#include "mrfbuildcore.h"

#include <getopt.h>

using std::cout;
using std::endl;

static const string option_message =
    "Input options:\n"
    " --msa <fmt>               choose format of MSA file\n"
    "                           fasta: aligned FASTA (default)\n"
    "                           a3m: HHsearch A3M format\n"
    " --edge <edge_file>        specify list of edges to determine MRF architecture\n"
    "\n"
    "Output options:\n"
    " -o <mrf_file>             write MRF model to file\n"
    "\n"
    "Preprocessing options:\n"
    " --seqwt <int>             sequence weighting\n"
    "                           0: no sequence weighting\n"
    "                           1: Henikoff's position-based weights (default)\n"
    " --effnum <int>            effective number of sequences\n"
    "                           0: no effective number (default)\n"
    "                           1: exponential of average entropy\n"
    "\n"
    "Regularization options:\n"
    " --regul <int>             regularization of node and edge weights\n"
    "                           0: no\n"
    "                           1: L2 regularization (default)\n"
    " --regnode-lambda <float>  weighting factor for node regularization (default: 0.01)\n"
    " --regedge-lambda <float>  weighting factor for edge regularization (default: 0.2)\n"
    " --regedge-scale <int>     scaling edge regularization term\n"
    "                           0: no\n"
    "                           1: yes (default)\n"
    "\n"
    "Optimization options:\n"
    " --delta <float>           minimum rate of decrease for objective function (default: 1e-4)\n"
    "\n"
    " -h, --help                show this help message\n";

MRFBuildCommandLine::MRFBuildCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    opt.out_filename = "";
    opt.build_opt.msa_fmt = AFASTA;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFBuildCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFBuildProcessor*)processor)->build(opt.msa_filename, opt.out_filename);
}

void MRFBuildCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " build <msa_file> [options]" << endl
         << endl
         << option_message
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
                if (parse_regul(optarg, opt.build_opt.parameterizer_opt.regul)) break;
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
                if (parse_msa_fmt(optarg, opt.build_opt.msa_fmt)) break;
                else return false;
            case 6:
                if (parse_str(optarg, opt.build_opt.eidx_filename)) break;
                else return false;
            case 7:
                if (parse_float(optarg, opt.build_opt.optim_opt.delta)) break;
                else return false;
            case 8:
                if (parse_int(optarg, opt.build_opt.msa_analyzer_opt.seq_wt)) break;
                else return false;
            case 9:
                if (parse_int(optarg, opt.build_opt.msa_analyzer_opt.eff_num)) break;
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
    if (opt.build_opt.parameterizer_opt.regul == RegulMethod::RegulMethod::L2) {
        opt.build_opt.parameterizer_opt.l2_opt.lambda1 = regnode_lambda;
        opt.build_opt.parameterizer_opt.l2_opt.lambda2 = regedge_lambda;
        opt.build_opt.parameterizer_opt.l2_opt.sc = regedge_scale;
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
