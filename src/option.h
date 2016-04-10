#ifndef _OPTION_H_
#define _OPTION_H_

#include <string>

#include "option.h"
#include "msaanalyze.h"
#include "parameterize.h"

using std::string;

enum SubCommand { NONE, HELP, BUILD, INFER, STAT, SHOW };

namespace Main {
    static const string command_list_message =
        "generate MRF model\n"
        "  build        Build MRF model\n"
        "\n"
        "examine evolutionary information\n"
        "  infer        Estimate sequence distribution with MRF model\n"
        "  stat         Estimate evolutionary constraints for MRF model\n"
        "  show         Show MRF model parameters\n"
        "\n"
        "Use `pmrf --version' to print the PMRF version\n"
        "Use `pmrf -h' or `pmrf --help' to print this help message\n"
        "Use `pmrf <command> -h' to print the help message on a specific command\n";
}

namespace Build {
    struct Option {
        string msa_filename;
        string out_filename;
        string eidx_filename;
        MSAFormat msa_fmt;
        MSAAnalyzer::Option msa_analyzer_opt;
        MRFParameterizer::ObjectiveFunction::Option parameterizer_opt;
        MRFParameterizer::Parameter::Option optim_opt;
    };

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
}

namespace Stat {
    enum Mode { MODE_PAIR, MODE_POS };

    enum Correct { CORR_NONE = 0, CORR_APC = 1, CORR_NCPS = 2 };

    class Option {
      public:
        Option(const Mode& mode=MODE_PAIR, const Correct& corr=CORR_APC, const bool& zscore=true)
        : mode(mode), corr(corr), zscore(zscore) {};

        string mrf_filename;
        Mode mode;
        Correct corr;
        bool zscore;
    };

    static const string option_message =
        "Options:\n"
        " --mode <mode>             calculation mode of evolutionary constraints\n"
        "                           pos: positional mode\n"
        "                           pair: pairwise mode (default)\n"
        " --corr <int>              pairwise coevolution score correction\n"
        "                           0: no correction\n"
        "                           1: average product (default)\n"
        "                           2: normalized coevolutionary pattern similarity\n"
        "\n"
        " -h, --help                show this help message\n";
}

namespace Infer {
    struct Option {
        string mrf_filename;
        string seq_filename;
    };

    static const string option_message =
        "Options:\n"
        " -h, --help                show this help message\n";
}

namespace Show {
    struct Option {
        string mrf_filename;
    };

    static const string option_message =
        "Options:\n"
        " -h, --help                show this help message\n";
}

#endif
