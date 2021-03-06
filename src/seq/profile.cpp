#include "profile.h"

#include "bgfreq.h"

static const NullSeqWeightEstimator NULL_SEQ_WEIGHT_ESTIMATOR;
static const NullEffSeqNumEstimator NULL_EFF_SEQ_NUM_ESTIMATOR;
static const NullMSAFilter NULL_MSA_FILTER;
static const MLEmitProbEstimator NULL_EMIT_PROB_ESTIMATOR;

Profile::Profile(const size_t& length, const Alphabet& abc) : length(length), abc(abc) {
    seq = "";
    char c = abc.get_unknown()[0];
    for (size_t i = 0; i < length; ++i) seq += c;
    init();
}

Profile::Profile(const string& seq, const Alphabet& abc) : seq(seq), abc(abc) {
    length = seq.size();
    init();
}

void Profile::init() {
    prob.resize(length, abc.get_canonical_size());
    eff_num.resize(length);
    bgfreq.resize(abc.get_canonical_size());
}

MatrixXf Profile::get_ll() const {
    int n = prob.rows();
    int m = prob.cols();
    MatrixXf s(n, m);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            s(i, j) = log(prob(i, j) / bgfreq(j));
    return s;
}

ProfileBuilder::ProfileBuilder(const Alphabet& abc) : abc(abc) {
    seq_weight_estimator = &NULL_SEQ_WEIGHT_ESTIMATOR;
    eff_seq_num_estimator = &NULL_EFF_SEQ_NUM_ESTIMATOR;
    msa_filter = &NULL_MSA_FILTER;
    emit_prob_estimator = &NULL_EMIT_PROB_ESTIMATOR;
}

Profile ProfileBuilder::build(const TraceVector& traces) const {
    int n = traces[0].get_length();
    Profile profile(n, abc);
    for (int i = 0; i < n; ++i) {
        TraceVector trs;
        vector<string> msa;
        collect_trace(traces, i, trs, msa);
        VectorXf wt = seq_weight_estimator->estimate(msa);
        double eff_num = eff_seq_num_estimator->estimate(msa, wt);
        VectorXf f = calc_freq(trs, i, wt, eff_num);    // frequency
        if (f.sum() > 0) profile.set_prob(i, emit_prob_estimator->estimate(f));
        else profile.set_prob(i, ROBINSON_BGFREQ.get_array(abc));
        profile.set_eff_num(i, eff_num);
        profile.set_bgfreq(ROBINSON_BGFREQ.get_array(abc));
    }
    return profile;
}

void ProfileBuilder::collect_trace(const TraceVector& traces, const size_t& idx, TraceVector& r_trs, vector<string>& r_msa) const {
    r_trs = traces.subset_matched(idx);
    r_msa = r_trs.get_matched_aseq_vec();
    if (r_msa.empty()) return;
    r_msa = msa_filter->filter(r_msa);
}

VectorXf ProfileBuilder::calc_freq(const TraceVector& trs, const size_t& idx, const VectorXf& wt, const double& eff_num) const {
    size_t n = trs.size();
    VectorXf v = VectorXf::Zero(abc.get_canonical_size());
    for (size_t i = 0; i < n; ++i)
        v += abc.get_count(trs[i].get_symbol_at(idx)) * wt(i);
    v /= v.sum();
    return v * eff_num;
}
