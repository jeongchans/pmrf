#ifndef _MSAANALYZE_H_
#define _MSAANALYZE_H_

#include <memory>

#include "seq/seqweight.h"
#include "seq/effseqnum.h"
#include "seq/msafilter.h"

using std::shared_ptr;

typedef int SeqWeightMethod;
const int NO_WEIGHT = 0;
const int POSITION_BASED = 1;

typedef int EffSeqNumMethod;
const int NO_EFFNUM = 0;
const int EXP_ENTROPY = 1;

class MSAAnalyzer {
  public:

    class Option {
      public:
        Option(const SeqWeightMethod& seq_wt=POSITION_BASED, const EffSeqNumMethod& eff_num=NO_EFFNUM)
        : seq_wt(seq_wt),
          eff_num(eff_num) {};

        int seq_wt;
        int eff_num;
    };

    MSAAnalyzer(Option& opt, const Alphabet& abc);

    shared_ptr<SeqWeightEstimator> seq_weight_estimator;
    shared_ptr<EffSeqNumEstimator> eff_seq_num_estimator;
    shared_ptr<TerminalGapRemover> termi_gap_remover;
};

#endif
