#include "msaanalyze.h"

MSAAnalyzer::MSAAnalyzer(Option& opt, const Alphabet& abc) {
    if (opt.seq_wt == MSAProcOption::SW_NO) {
        seq_weight_estimator = std::shared_ptr<NullSeqWeightEstimator>(
            new NullSeqWeightEstimator());
    } else if (opt.seq_wt == MSAProcOption::SW_PB) {
        seq_weight_estimator = std::shared_ptr<PBSeqWeightEstimator>(
            new PBSeqWeightEstimator());
    } else if (opt.seq_wt == MSAProcOption::SW_CLSTR) {
        seq_weight_estimator = std::shared_ptr<ClstrSeqWeightEstimator>(
            new ClstrSeqWeightEstimator(abc, opt.clstr_maxidt));
    }

    if (opt.eff_num == MSAProcOption::NEFF_NO) {
        eff_seq_num_estimator = std::shared_ptr<NullEffSeqNumEstimator>(
            new NullEffSeqNumEstimator());
    } else if (opt.eff_num == MSAProcOption::NEFF_SHANNON) {
        eff_seq_num_estimator = std::shared_ptr<ExpEntropyEffSeqNumEstimator>(
            new ExpEntropyEffSeqNumEstimator(abc, seq_weight_estimator.get()));
    } else if (opt.eff_num == MSAProcOption::NEFF_CLSTR) {
        eff_seq_num_estimator = std::shared_ptr<ClstrEffSeqNumEstimator>(
            new ClstrEffSeqNumEstimator(abc, opt.clstr_maxidt));
    }

    termi_gap_remover = std::shared_ptr<TerminalGapRemover>(
        new TerminalGapRemover(abc, opt.termi_maxgapperc));
}
