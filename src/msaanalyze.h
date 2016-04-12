#ifndef _MSAANALYZE_H_
#define _MSAANALYZE_H_

#include <memory>

#include "seq/seqweight.h"
#include "seq/effseqnum.h"
#include "seq/msafilter.h"

using std::shared_ptr;

namespace MSAProcOption {
    enum SeqWeight { SW_NO = 0, SW_PB = 1 };

    enum EffSeqNum { NEFF_NO = 0, NEFF_ENTROPY = 1, NEFF_CLSTR = 2 };
}

class MSAAnalyzer {
  public:

    class Option {
      public:
        Option(const MSAProcOption::SeqWeight& seq_wt=MSAProcOption::SW_PB, const MSAProcOption::EffSeqNum& eff_num=MSAProcOption::NEFF_CLSTR)
        : seq_wt(seq_wt),
          eff_num(eff_num) {};

        MSAProcOption::SeqWeight seq_wt;
        MSAProcOption::EffSeqNum eff_num;
    };

    MSAAnalyzer(Option& opt, const Alphabet& abc);

    shared_ptr<SeqWeightEstimator> seq_weight_estimator;
    shared_ptr<EffSeqNumEstimator> eff_seq_num_estimator;
    shared_ptr<TerminalGapRemover> termi_gap_remover;
};

#endif
