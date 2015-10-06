#include "mrfbuildcommand.h"
#include "mrfbuildcore.h"

#include <getopt.h>

using std::cout;
using std::endl;

static const string option_message =
    "Input options:\n"
    " --msa <fmt>               Multiple sequence alignment format\n"
    "                           fasta: aligned FASTA (default)\n"
    "                           a3m: HHsearch A3M format\n"
    " --edge <file>             Restrict edges to list of index pairs\n"
    "\n"
    "Preprocessing options:\n"
    " --seqwt <int>             Sequence weight method\n"
    "                           0: no sequence weight\n"
    "                           1: Henikoff's position-based weights (default)\n"
    "\n"
    "Regularization options:\n"
    " --regnode <int>           Regularization of node weights\n"
    "                           0: no\n"
    "                           1: L2 regularization (default)\n"
    "                           2: profile-based regularization\n"
    " --regnode-lambda <float>  Weighting factor for node regularization (default: 0.01)\n"
    "\n"
    " --regedge <int>           Regularization of edge weights\n"
    "                           0: none\n"
    "                           1: L2 regularization (default)\n"
    " --regedge-lambda <float>  Weighting factor for edge regularization (default: 0.2)\n"
    " --regedge-scale <int>     Scaling edge regularization\n"
    "                           0: no\n"
    "                           1: yes (default)\n"
    "\n"
    "Optimization options:\n"
    " --delta <float>           Minimum rate of decrease of objective function (default: 1e-4)\n"
    "\n"
    " -h, --help                Help\n";

MRFBuildCommandLine::MRFBuildCommandLine(int argc, char** argv) : MRFCommandLine(argc, argv) {
    opt.out_filename = "";
    opt.build_opt.msa_fmt = AFASTA;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFBuildCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFBuildProcessor*)processor)->build(opt.msa_filename, opt.out_filename);
}

void MRFBuildCommandLine::show_help() {
    cout << "Usage: " << PROGNAME << " build <file> [options]" << endl
         << endl
         << option_message
         << endl;
}

bool MRFBuildCommandLine::parse_command_line(int argc, char** argv) {
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help", 0, 0, 0},
        {"regnode", required_argument, 0, 0},
        {"regnode-lambda", required_argument, 0, 0},
        {"regedge", required_argument, 0, 0},
        {"regedge-lambda", required_argument, 0, 0},
        {"regedge-scale", required_argument, 0, 0},
        {"msa", required_argument, 0, 0},
        {"edge", required_argument, 0, 0},
        {"delta", required_argument, 0, 0},
        {"seqwt", required_argument, 0, 0},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    double node_regul_lambda, edge_regul_lambda;
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
                if (parse_regnode(optarg, opt.build_opt.parameterizer_opt.node_regul)) break;
                else return false;
            case 2:
                if (parse_double(optarg, node_regul_lambda)) break;
                else return false;
            case 3:
                if (parse_regedge(optarg, opt.build_opt.parameterizer_opt.edge_regul)) break;
                else return false;
            case 4:
                if (parse_double(optarg, edge_regul_lambda)) break;
                else return false;
            case 5:
                if (parse_bool(optarg, opt.build_opt.parameterizer_opt.edge_l2_opt.sc)) break;
                else return false;
            case 6:
                if (parse_msa_fmt(optarg, opt.build_opt.msa_fmt)) break;
                else return false;
            case 7:
                if (parse_str(optarg, opt.build_opt.eidx_filename)) break;
                else return false;
            case 8:
                if (parse_float(optarg, opt.build_opt.optim_opt.delta)) break;
                else return false;
            case 9:
                if (parse_int(optarg, opt.build_opt.msa_analyzer_opt.seq_wt)) break;
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
            error_message = "Unknown option\n";
            return false;
        }
    }
    optind++;
    if (optind < argc) {
        opt.msa_filename = argv[optind];
    } else {
        error_message = "Not enough arguments\n";
        return false;
    }
    if (opt.build_opt.parameterizer_opt.node_regul == NodeRegulMethod::L2)
        opt.build_opt.parameterizer_opt.node_l2_opt.lambda = node_regul_lambda;
    else if (opt.build_opt.parameterizer_opt.node_regul == NodeRegulMethod::PSSM)
        opt.build_opt.parameterizer_opt.node_pssm_opt.lambda = node_regul_lambda;
    if (opt.build_opt.parameterizer_opt.edge_regul == EdgeRegulMethod::L2)
        opt.build_opt.parameterizer_opt.edge_l2_opt.lambda = edge_regul_lambda;
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

bool MRFBuildCommandLine::parse_regnode(char* optarg, NodeRegulMethod& arg) {
    int d = atoi(optarg);
    if (d == 0) arg = NodeRegulMethod::NONE;
    else if (d == 1) arg = NodeRegulMethod::L2;
    else if (d == 2) arg = NodeRegulMethod::PSSM;
    else {
        error_message = "Unknown node regularization option: " + d;
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::parse_double(char* optarg, double& arg) {
    arg = atof(optarg);
    return true;
}

bool MRFBuildCommandLine::parse_regedge(char* optarg, EdgeRegulMethod& arg) {
    int d = atoi(optarg);
    if (d == 0) arg = EdgeRegulMethod::NONE;
    else if (d == 1) arg = EdgeRegulMethod::L2;
    else {
        error_message = "Unknown edge regularization option: " + d;
        return false;
    }
    return true;
}
