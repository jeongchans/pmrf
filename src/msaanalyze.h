#ifndef _MSAANALYZE_H_
#define _MSAANALYZE_H_

#include <memory>

#include "seq/seqweight.h"
#include "seq/effseqnum.h"
#include "seq/msafilter.h"
#include "seq/trace.h"

using std::shared_ptr;

namespace MSAProcOption {
    enum SeqWeight { SW_NO = 0, SW_PB = 1, SW_CLSTR = 2 };

    enum EffSeqNum { NEFF_NO = 0, NEFF_SHANNON = 1, NEFF_CLSTR = 2, NEFF_JOINT_SHANNON = 3 };
}

class MSAAnalyzer {
  public:

    class Option {
      public:
        Option()
        : seq_wt(MSAProcOption::SW_CLSTR),
          eff_num(MSAProcOption::NEFF_CLSTR),
          clstr_maxidt(0.8),
          termi_maxgapperc(0.1) {};

        MSAProcOption::SeqWeight seq_wt;
        MSAProcOption::EffSeqNum eff_num;

        float clstr_maxidt;
        float termi_maxgapperc;
    };

    MSAAnalyzer(const Option& opt, const Alphabet& abc);

    void calc_sw_and_neff(const TraceVector& traces, VectorXf& sw, float& neff) const;

    shared_ptr<SeqWeightEstimator> seq_weight_estimator;
    shared_ptr<EffSeqNumEstimator> eff_seq_num_estimator;
    shared_ptr<TerminalGapRemover> termi_gap_remover;

  private:
    const Option& opt;
};

#endif
