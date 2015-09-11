#ifndef _MSAANALYZE_H_
#define _MSAANALYZE_H_

#include <memory>

#include "msaanalysis/seqweight.h"
#include "msaanalysis/effseqnum.h"
#include "msaanalysis/msafilter.h"

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

    std::shared_ptr<SeqWeightEstimator> seq_weight_estimator;
    std::shared_ptr<EffSeqNumEstimator> eff_seq_num_estimator;
    std::shared_ptr<TerminalGapRemover> termi_gap_remover;

  private:
    const Alphabet& abc;
};

#endif
