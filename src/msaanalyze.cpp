#include "msaanalyze.h"

using std::make_shared;

MSAAnalyzer::MSAAnalyzer(const Option& opt, const Alphabet& abc) : opt(opt) {
    if (opt.seq_wt == MSAProcOption::SW_NO) {
        seq_weight_estimator = make_shared<NullSeqWeightEstimator>();
    } else if (opt.seq_wt == MSAProcOption::SW_PB) {
        seq_weight_estimator = make_shared<PBSeqWeightEstimator>();
    } else if (opt.seq_wt == MSAProcOption::SW_CLSTR) {
        seq_weight_estimator = make_shared<ClstrSeqWeightEstimator>(abc, opt.clstr_maxidt);
    }

    if (opt.eff_num == MSAProcOption::NEFF_NO) {
        eff_seq_num_estimator = make_shared<NullEffSeqNumEstimator>();
    } else if (opt.eff_num == MSAProcOption::NEFF_SHANNON) {
        eff_seq_num_estimator = make_shared<ExpEntropyEffSeqNumEstimator>(abc, seq_weight_estimator.get());
    } else if (opt.eff_num == MSAProcOption::NEFF_JOINT_SHANNON) {
        eff_seq_num_estimator = make_shared<ExpJointEntropyEffSeqNumEstimator>(abc, seq_weight_estimator.get());
    } else if (opt.eff_num == MSAProcOption::NEFF_CLSTR) {
        if (opt.seq_wt == MSAProcOption::SW_CLSTR) {
            eff_seq_num_estimator = make_shared<NullEffSeqNumEstimator>();
        } else {
            eff_seq_num_estimator = make_shared<ClstrEffSeqNumEstimator>(abc, opt.clstr_maxidt);
        }
    }

    termi_gap_remover = make_shared<TerminalGapRemover>(abc, opt.termi_maxgapperc);
}

void MSAAnalyzer::calc_sw_and_neff(const TraceVector& traces, VectorXf& sw, float& neff) const {
    vector<string> msa = traces.get_matched_aseq_vec();
    msa = termi_gap_remover->filter(msa);
    sw = seq_weight_estimator->estimate(msa);
    neff = 1.;  // neff ignored
    //TODO
    //if (opt.seq_wt == MSAProcOption::SW_CLSTR && opt.eff_num == MSAProcOption::NEFF_CLSTR) neff = sw.sum();
    //else neff = eff_seq_num_estimator->estimate(msa, sw / sw.sum());
}
