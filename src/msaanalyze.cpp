#include "msaanalyze.h"

MSAAnalyzer::MSAAnalyzer(Option& opt, const Alphabet& abc) {
    if (opt.seq_wt == NO_WEIGHT) {
        seq_weight_estimator = std::shared_ptr<NullSeqWeightEstimator>(
            new NullSeqWeightEstimator());
    } else if (opt.seq_wt == POSITION_BASED) {
        seq_weight_estimator = std::shared_ptr<PBSeqWeightEstimator>(
            new PBSeqWeightEstimator());
    }
    if (opt.eff_num == NO_EFFNUM) {
        eff_seq_num_estimator = std::shared_ptr<NullEffSeqNumEstimator>(
            new NullEffSeqNumEstimator());
    } else if (opt.eff_num == EXP_ENTROPY) {
        eff_seq_num_estimator = std::shared_ptr<ExpEntropyEffSeqNumEstimator>(
            new ExpEntropyEffSeqNumEstimator(abc, seq_weight_estimator.get()));
    }
    termi_gap_remover = std::shared_ptr<TerminalGapRemover>(
        new TerminalGapRemover(abc, 0.1));
}
