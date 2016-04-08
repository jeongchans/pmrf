#ifndef _OPTION_H_
#define _OPTION_H_

#include <string>

#include "option.h"
#include "msaanalyze.h"
#include "parameterize.h"

using std::string;

enum SubCommand { NONE, HELP, BUILD, INFER, STAT };

namespace Build {
    struct Option {
        MSAFormat msa_fmt;
        string eidx_filename;
        MSAAnalyzer::Option msa_analyzer_opt;
        MRFParameterizer::ObjectiveFunction::Option parameterizer_opt;
        MRFParameterizer::Parameter::Option optim_opt;
    };
}

namespace Stat {
    enum Mode { MODE_PAIR, MODE_POS };

    enum Correct { CORR_NONE = 0, CORR_APC = 1, CORR_NCPS = 2 };

    class Option {
      public:
        Option(const Mode& mode=MODE_PAIR, const Correct& corr=CORR_APC, const bool& zscore=true)
        : mode(mode), corr(corr), zscore(zscore) {};

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
}

#endif
