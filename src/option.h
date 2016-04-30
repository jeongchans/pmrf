#ifndef _OPTION_H_
#define _OPTION_H_

#include <string>

#include "msaanalyze.h"
#include "parameterize.h"

using std::string;

enum SubCommand { NONE, HELP, BUILD, INFER, STAT, SHOW };

namespace Main {
    static const string usage_message = "<command> [<args>]";
    static const string command_list_message =
        "The following commands are available for MRF modeling and varisous applications:\n"
        "\n"
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
        MRFParameterizer::Option parameterizer_opt;
        MRFParameterizer::Parameter::Option optim_opt;
    };

    static const string usage_message = "build <msa_file> [options]";
    static const string option_message =
        "Input options:\n"
        " --msa <fmt>               choose format of MSA file\n"
        "                           fasta: aligned FASTA (default)\n"
//        "                           flat: list of aligned sequences\n"
        "                           a3m: HHsearch A3M format\n"
//        " --ref <int>               index of reference sequence (default: 1)\n"
//        "                           0 or less number will consider all the MSA columns\n"
        " -e, --edge <edge_file>    specify list of edges to determine MRF architecture\n"
        "\n"
        "Output options:\n"
        " -o, --output <mrf_file>   write MRF model to file\n"
        "\n"
        "Preprocessing options:\n"
        " --seqwt <method>          sequence weighting\n"
        "                           no: no sequence weighting\n"
//        "                           pb: Henikoff's position-based weights (default)\n"
        "                           clstr: Henikoff's position-based weights (default)\n"
//        " --neff <method>           effective number of sequences\n"
//        "                           no: no effective number\n"
//        "                           clstr: number of sequences clustered by sequence identity (default)\n"
//        "                           shan: exponential of average shannon entropy of amino acid pair\n"
        " --clstr-maxidt <float>    maximum sequence identity between sequence clusters (default: 0.8)\n"
        "\n"
        "Parameterization options:\n"
        " --symmetric               use symmetry constraint on edge weights\n"
        " --regul <method>          regularization of node and edge weights\n"
        "                           no: no\n"
        "                           l2: L2 regularization (default)\n"
        " --regul-vl <float>        weighting factor for node regularization (default: auto)\n"
        " --regul-wl <float>        weighting factor for edge regularization (default: auto)\n"
        //" --regw-sc-deg <yes|no>    scaling edge regularization w.r.t. average degree (default: yes)\n"
        //" --regw-sc-neff <yes|no>   scaling edge regularization w.r.t. effective number of sequences\n"
        //"                           (default: yes)\n"
//        " --no-regw-sc-deg          disable edge regularization scaling w.r.t. average degree\n"
//        " --no-regw-sc-neff         disable edge regularization scaling w.r.t. effective number of sequences\n"
        "\n"
        "Optimization options:\n"
        " --lbfgs-corr <int>        number of correction (default: 100)\n"
        " --lbfgs-epsilon <float>   convergence criterion (default: 1e-5)\n"
        " --lbfgs-delta <float>     stopping criterion (default: 1e-7)\n"
        " --lbfgs-maxiter <int>     maximum iteration (default: 500)\n"
        "\n"
        " -h, --help                show this help message\n";
}

namespace Stat {
    enum Mode { MODE_PAIR, MODE_POS };
    enum Correct { CORR_NONE, CORR_APC, CORR_NCPS };

    class Option {
      public:
        Option(const Mode& mode=MODE_PAIR, const Correct& corr=CORR_APC, const bool& zscore=true)
        : mode(mode), corr(corr), zscore(zscore) {};

        string mrf_filename;
        Mode mode;
        Correct corr;
        bool zscore;
    };

    static const string usage_message = "stat <mrf_file> [options]";
    static const string option_message =
        "Options:\n"
        " --mode <mode>             calculation mode of evolutionary constraints\n"
        "                           pos: positional mode\n"
        "                           pair: pairwise mode (default)\n"
        "\n"
        "Pairwise coevolution options:\n"
        " --corr <method>           pairwise coevolution score correction\n"
        "                           no: no correction\n"
        "                           apc: average product (default)\n"
        "                           ncps: normalized coevolutionary pattern similarity\n"
        "\n"
        " -h, --help                show this help message\n";
}

namespace Infer {
    struct Option {
        string mrf_filename;
        string seq_filename;
    };

    static const string usage_message = "infer <mrf_file> <seq_file> [options]";
    static const string option_message =
        "Options:\n"
//        " --ref <seq_file>          calculate the score difference for <seq_file> (FASTA)\n"
//        " --offset <float>          score offset (default: 0.0)\n"
        "\n"
        " -h, --help                show this help message\n";
}

namespace Show {
    struct Option {
        string mrf_filename;
    };

    static const string usage_message = "show <mrf_file> [options]";
    static const string option_message =
        "Options:\n"
//        " --brief                   show only the brief information\n"
//        " --list-prof               list profile parameters\n"
//        " --list-mrf                list MRF parameters\n"
        "\n"
        " -h, --help                show this help message\n";
}

#endif
