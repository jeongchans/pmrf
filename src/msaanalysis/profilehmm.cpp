#include "profilehmm.h"

static const NullSeqWeightEstimator NULL_SEQ_WEIGHT_ESTIMATOR;
static const NullEffSeqNumEstimator NULL_EFF_SEQ_NUM_ESTIMATOR;
static const NullMSAFilter NULL_MSA_FILTER;
static const MLEmitProbEstimator NULL_EMIT_PROB_ESTIMATOR;
static const MLTransitProbEstimator NULL_TRANSIT_PROB_ESTIMATOR;

ProfileHMM::ProfileHMM(const size_t& length, const Alphabet& abc) : length(length), abc(abc) {
    seq = "";
    char c = abc.get_unknown()[0];
    for (size_t i = 0; i < length; ++i) seq += c;
    init();
}

ProfileHMM::ProfileHMM(const std::string& seq, const Alphabet& abc) : seq(seq), abc(abc) {
    length = seq.size();
    init();
}

void ProfileHMM::init() {
    m_state.clear();
    d_state.clear();
    i_state.clear();
    for (size_t i = 0; i < length; ++i) {
        m_state.push_back(HMMMatchState(abc.get_canonical_size()));
        d_state.push_back(HMMDeleteState());
        i_state.push_back(HMMInsertState(abc.get_canonical_size()));
    }
    null_emit.resize(abc.get_canonical_size());
    null_emit = 0;
    for (size_t i = 0; i < NUM_STATE_TYPE; ++i) {
        null_transit[i].resize(NUM_STATE_TYPE);
        null_transit[i] = 0;
    }
}

void ProfileHMM::set_seq(const std::string& seq) {
    this->seq = seq;
    length = seq.size();
}

void ProfileHMM::traverse(HMMStateVisitor& visitor) {
    for (size_t i = 0; i < length; ++i) {
        m_state[i].accept(&visitor, i);
        d_state[i].accept(&visitor, i);
        i_state[i].accept(&visitor, i);
    }
}

HMMParamEstimator::HMMParamEstimator(const TraceVector& traces, const Alphabet& abc, const bool& termi_gap_filter) : traces(traces), abc(abc), termi_gap_filter(termi_gap_filter) {
    seq_weight_estimator = &NULL_SEQ_WEIGHT_ESTIMATOR;
    eff_seq_num_estimator = &NULL_EFF_SEQ_NUM_ESTIMATOR;
    emit_prob_estimator = &NULL_EMIT_PROB_ESTIMATOR;
    transit_prob_estimator = &NULL_TRANSIT_PROB_ESTIMATOR;
    msa_filter = &NULL_MSA_FILTER;
}

inline Float1dArray HMMParamEstimator::calc_seq_weight(const std::vector<std::string>& msa) {
    return seq_weight_estimator->estimate(msa);
}

inline double HMMParamEstimator::calc_eff_seq_num(const std::vector<std::string>& msa) {
    return eff_seq_num_estimator->estimate(msa);
}

inline Float1dArray HMMParamEstimator::calc_emit_freq(const StateType& type, const size_t& idx, const TraceVector& trs, const Float1dArray& wt, const double& eff_num) {
    size_t n = trs.size();
    Float1dArray v(abc.get_canonical_size());
    v = 0;
    std::string s;
    for (size_t i = 0; i < n; ++i) {
        s = trs[i].get_emit(type, idx);
        for (std::string::iterator pos = s.begin(); pos != s.end(); ++pos)
            v += abc.get_count(*pos) * wt(i);
    }
    scale(v, eff_num);
    return v;
}

inline Float1dArray HMMParamEstimator::calc_tr_freq(const StateType& type, const size_t& idx, const TraceVector& trs, const Float1dArray& wt, const double& eff_num) {
    size_t n = trs.size();
    Float1dArray v(NUM_STATE_TYPE);
    v = 0;
    for (size_t i = 0; i < n; ++i)
        v += trs[i].get_transit_count(type, idx) * wt(i);
    scale(v, eff_num);
    return v;
}

inline Float1dArray HMMParamEstimator::estimate_emit_prob(const Float1dArray& f) {
    return emit_prob_estimator->estimate(f);
}

inline Float1dArray HMMParamEstimator::estimate_transit_prob(const Float1dArray& f, const StateType& type) {
    return transit_prob_estimator->estimate(f, type);
}

inline void HMMParamEstimator::collect_trace(const StateType type, const size_t& idx, TraceVector& trs, std::vector<std::string>& msa) {
    trs = traces.subset_passing(type, idx, termi_gap_filter);
    msa = trs.get_MD_seq_vec();
    if (msa.empty()) return;
    msa = msa_filter->filter(msa);
}

void HMMParamEstimator::visit_match(HMMMatchState* state, const size_t& idx) {
    TraceVector m_traces;
    std::vector<std::string> msa;
    collect_trace(MATCH, idx, m_traces, msa);
    // estimate frequencies
    Float1dArray f_v = zeros(abc.get_canonical_size());
    Float1dArray t_v = zeros(NUM_STATE_TYPE);
    double eff_num = 0;
    if (!m_traces.empty()) {
        Float1dArray wt = calc_seq_weight(msa);
        eff_num = calc_eff_seq_num(msa);
        f_v = calc_emit_freq(MATCH, idx, m_traces, wt, eff_num); // emit frequency
        t_v = calc_tr_freq(MATCH, idx, m_traces, wt, eff_num);   // transit frequency
    }
    // estimate probablities from frequencies and set HMM parameters
    state->set_emit(estimate_emit_prob(f_v));
    state->set_transit(estimate_transit_prob(t_v, MATCH));
    state->set_eff_num(eff_num);
}

void HMMParamEstimator::visit_delete(HMMDeleteState* state, const size_t& idx) {
    TraceVector m_traces;
    std::vector<std::string> msa;
    collect_trace(DELETE, idx, m_traces, msa);
    // estimate frequencies
    Float1dArray t_v = zeros(NUM_STATE_TYPE);
    double eff_num = 0;
    if (!m_traces.empty()) {
        Float1dArray wt = calc_seq_weight(msa);
        eff_num = calc_eff_seq_num(msa);
        t_v = calc_tr_freq(DELETE, idx, m_traces, wt, eff_num);
    }
    // estimate probablities from frequencies and set HMM parameters
    state->set_transit(estimate_transit_prob(t_v, DELETE));
    state->set_eff_num(eff_num);
}

void HMMParamEstimator::visit_insert(HMMInsertState* state, const size_t& idx) {
    TraceVector m_traces;
    std::vector<std::string> msa;
    collect_trace(INSERT, idx, m_traces, msa);
    // estimate frequencies
    Float1dArray f_v = zeros(abc.get_canonical_size());
    Float1dArray t_v = zeros(NUM_STATE_TYPE);
    double eff_num = 0;
    if (!m_traces.empty()) {
        Float1dArray wt = calc_seq_weight(msa);
        eff_num = calc_eff_seq_num(msa);
        f_v = calc_emit_freq(INSERT, idx, m_traces, wt, eff_num);
        t_v = calc_tr_freq(INSERT, idx, m_traces, wt, eff_num);
    }
    // estimate probablities from frequencies and set HMM parameters
    state->set_emit(estimate_emit_prob(f_v));
    state->set_transit(estimate_transit_prob(t_v, INSERT));
    state->set_eff_num(eff_num);
}

void HMMParamEstimator::parameterize(ProfileHMM& model) {
    vector<string> msa = traces.get_MD_seq_vec();
    msa = msa_filter->filter(msa);
    model.set_eff_num(eff_seq_num_estimator->estimate(msa));
    model.traverse(*this);
}

void PseudocountMixer::visit_match(HMMMatchState* state, const size_t&) {
    double eff_num = state->get_eff_num();
    Float1dArray em = state->get_emit().copy();
    em *= eff_num;
    state->set_emit(emit_prob_estimator->estimate(em));
    Float1dArray tr = state->get_transit().copy();
    tr *= eff_num;
    state->set_transit(transit_prob_estimator->estimate(tr, MATCH));
}

void PseudocountMixer::visit_delete(HMMDeleteState* state, const size_t&) {
    double eff_num = state->get_eff_num();
    Float1dArray tr = state->get_transit().copy();
    tr *= eff_num;
    state->set_transit(transit_prob_estimator->estimate(tr, DELETE));
}

void PseudocountMixer::visit_insert(HMMInsertState* state, const size_t&) {
    double eff_num = state->get_eff_num();
    Float1dArray tr = state->get_transit().copy();
    tr *= eff_num;
    state->set_transit(transit_prob_estimator->estimate(tr, INSERT));
}
