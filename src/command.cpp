#include "command.h"

#include <getopt.h>

#include "core.h"

using std::cout;
using std::endl;

MRFCommandLine::MRFCommandLine(int argc, char** argv, const string& usage_message, const string& opt_message) 
: argc(argc), 
  usage_message(usage_message), 
  opt_message(opt_message),
  validity(false), 
  error_message("") {
    this->argv = (char**)malloc(sizeof(char*) * argc);
    for (int i = 0; i < argc; ++i) {
        this->argv[i] = (char*)malloc(strlen(argv[i]) + 1);
        strcpy(this->argv[i], argv[i]);
    }
}

MRFCommandLine::~MRFCommandLine() {
    for (int i = 0; i < argc; ++i) free(argv[i]);
    free(argv);
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

MRFMainCommandLine::MRFMainCommandLine(int argc, char** argv) 
: MRFCommandLine(argc, argv, Main::usage_message, Main::command_list_message) {
    subcmd = NONE;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFMainCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFMainProcessor*)processor)->run_mrf_cmd(subcmd);
}

bool MRFMainCommandLine::parse_command_line(int argc, char** argv) {
    if (argc < 2) { show_help(); exit(0); }
    if (argv[1][0] != '-') {
        string cmd = argv[1];
        if (cmd == "help") subcmd = HELP;
        else if (cmd == "build") subcmd = BUILD;
        else if (cmd == "infer") subcmd = INFER;
        else if (cmd == "stat") subcmd = STAT;
        else if (cmd == "show") subcmd = SHOW;
        else {
            error_message = string("Unknown command: '") + cmd + "'\n" + "Use `pmrf --help'";
            return false;
        }
        return true;
    }
    static struct option opts[] = {
        {"help",            no_argument,        0, 'h'},
        {"version",         no_argument,        0, 100},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        else if (c == 100) { show_version(); exit(0); }
        else if (c == 'h') { show_help(); exit(0); }
        else return false;
    }
    return true;
}

/**
   @class MRFBuildCommandLine
 */

MRFBuildCommandLine::MRFBuildCommandLine(int argc, char** argv) 
: MRFCommandLine(argc, argv, Build::usage_message, Build::option_message) {
    opt.out_filename = "";
    opt.msa_fmt = AFASTA;
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFBuildCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFBuildProcessor*)processor)->build();
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
        {"neff",            required_argument,  0, 301},
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
        //{"reg-lambda-c1",   required_argument,  0, -200},
        //{"reg-lambda-c2",   required_argument,  0, -201},
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
        else if (c == 301) validity = parse_str(optarg, sval) && set_eff_num(opt.msa_analyzer_opt.eff_num, sval);
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
        //else if (c == -200) validity = parse_float(optarg, opt.parameterizer_opt.reg_lambda_c1);
        //else if (c == -201) validity = parse_float(optarg, opt.parameterizer_opt.reg_lambda_c2);
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

bool MRFBuildCommandLine::set_eff_num(MSAProcOption::EffSeqNum& arg, const string& val) {
    if (val == "no") arg = MSAProcOption::NEFF_NO;
    else if (val == "clstr") arg = MSAProcOption::NEFF_CLSTR;
    else if (val == "shan2") arg = MSAProcOption::NEFF_JOINT_SHANNON;
    else {
        error_message = "Unsupported option for effective number estimation: " + val;
        return false;
    }
    return true;
}

bool MRFBuildCommandLine::set_regul(RegulMethod::RegulMethod& arg, const string& val) {
    if (val == "no") arg = RegulMethod::REGUL_NONE;
    else if (val == "l2") arg = RegulMethod::REGUL_L2;
    else return false;
    return true;
}

/**
   @class MRFStatCommandLine
 */

MRFStatCommandLine::MRFStatCommandLine(int argc, char** argv) 
: MRFCommandLine(argc, argv, Stat::usage_message, Stat::option_message) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFStatCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFStatProcessor*)processor)->stat();
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

MRFInferCommandLine::MRFInferCommandLine(int argc, char** argv) 
: MRFCommandLine(argc, argv, Infer::usage_message, Infer::option_message) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFInferCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFInferProcessor*)processor)->infer();
}

bool MRFInferCommandLine::parse_command_line(int argc, char** argv) {
    float fval;
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help",            no_argument,        0, 'h'},
        {"node-offset",     required_argument,  0, 200},
        {"edge-offset",     required_argument,  0, 201},
        {"prof-offset",     required_argument,  0, 300},
        {"gap-score",       required_argument,  0, 301},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        else if (c == 200) validity = parse_float(optarg, fval) && set_val<float>(opt.node_offset, fval);
        else if (c == 201) validity = parse_float(optarg, fval) && set_val<float>(opt.edge_offset, fval);
        else if (c == 300) validity = parse_float(optarg, fval) && set_val<float>(opt.prof_offset, fval);
        else if (c == 301) validity = parse_float(optarg, fval) && set_val<float>(opt.gap_score, fval);
        else if (c == 'h') { show_help(); exit(0); }
        else validity = false;
        if (!validity) return set_opt_err_msg("--" + string(opts[opt_idx].name), optarg);
    }
    if (optind != argc - 2) { show_help(); exit(EXIT_FAILURE); }
    opt.mrf_filename = argv[optind++];
    opt.seq_filename = argv[optind++];
    return true;
}

/**
   @class MRFShowCommandLine
 */

MRFShowCommandLine::MRFShowCommandLine(int argc, char** argv) 
: MRFCommandLine(argc, argv, Show::usage_message, Show::option_message) {
    validity = parse_command_line(this->argc, (char**)(this->argv));
}

int MRFShowCommandLine::process_command(MRFCmdProcessor *processor) {
    return ((MRFShowProcessor*)processor)->show();
}

bool MRFShowCommandLine::parse_command_line(int argc, char** argv) {
    string sval;
    int ival;
    optind = 0;     // initialize getopt_long()
    static struct option opts[] = {
        {"help",            no_argument,        0, 'h'},
        {"seq",             no_argument,        0, 100},
        {"profile",         no_argument,        0, 101},
        {"mrf",             no_argument,        0, 102},
        {"node",            required_argument,  0, 200},
        {"edge",            required_argument,  0, 201},
        {0, 0, 0, 0}
    };
    int opt_idx = 0;
    int c;
    while (true) {
        c = getopt_long(argc, argv, "h", opts, &opt_idx);
        if (c == -1) break;
        else if (c == 100) validity = set_val<bool>(opt.seq_flag, true);
        else if (c == 101) validity = set_val<bool>(opt.prof_flag, true);
        else if (c == 102) validity = set_val<bool>(opt.mrf_flag, true);
        else if (c == 200) validity = parse_int(optarg, ival) && ival > 0 && set_val<size_t>(opt.v_pos, (size_t) ival) && set_val<bool>(opt.node_flag, true);
        else if (c == 201) validity = parse_str(optarg, sval) && set_edge_idx(opt.w_pos, sval) && set_val<bool>(opt.edge_flag, true);
        else if (c == 'h') { show_help(); exit(0); }
        else validity = false;
        if (!validity) return set_opt_err_msg("--" + string(opts[opt_idx].name), optarg);
    }
    if (optind != argc - 1) { show_help(); exit(EXIT_FAILURE); }
    opt.mrf_filename = argv[optind++];
    return true;
}

bool MRFShowCommandLine::set_edge_idx(EdgeIndex& arg, const string& val) {
    size_t pos = val.find(",");
    if (pos == string::npos) return false;
    arg.idx1 = std::stoi(val.substr(0, pos));
    arg.idx2 = std::stoi(val.substr(pos + 1));
    return true;
}
