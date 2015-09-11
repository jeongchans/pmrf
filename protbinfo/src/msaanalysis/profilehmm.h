#ifndef _PROFILEHMM_H_
#define _PROFILEHMM_H_

#include <vector>
#include <map>

#include "common/common.h"
#include "common/numeric.h"
#include "seq/alphabet.h"
#include "align/trace.h"
#include "profile.h"
#include "transitprob.h"
#include "common/heteromap.h"
#include "hmmstate.h"

class ProfileHMM {
  public:
    ProfileHMM(const size_t& length, const Alphabet& abc);
    ProfileHMM(const std::string& seq, const Alphabet& abc);

    // state getters
    HMMMatchState& get_match(const size_t& idx) { return m_state[idx]; }
    HMMDeleteState& get_delete(const size_t& idx) { return d_state[idx]; }
    HMMInsertState& get_insert(const size_t& idx) { return i_state[idx]; }
    const HMMMatchState& get_match(const size_t& idx) const
        { return m_state[idx]; }
    const HMMDeleteState& get_delete(const size_t& idx) const
        { return d_state[idx]; }
    const HMMInsertState& get_insert(const size_t& idx) const
        { return i_state[idx]; }

    // length, sequence, description, and eff_num
    size_t get_length() const { return length; }
    std::string get_seq() const { return seq; }
    void set_seq(const std::string& seq);
    std::string get_emit_symbol() const { return abc.get_canonical(); }
    std::string get_desc() const { return desc; }
    void set_desc(const std::string& desc) { this->desc = desc; }
    double get_eff_num() const { return eff_num; }
    void set_eff_num(const double& eff_num) { this->eff_num = eff_num; }

    // null emission and transition probabilities
    Float1dArray get_null_emit() const { return null_emit.copy(); }
    Float1dArray get_null_transit(const StateType& type) const
        { return null_transit[type].copy(); }
    void set_null_emit(const Float1dArray& v) { null_emit = v; }
    void set_null_transit(const StateType& type, const Float1dArray& v)
        { null_transit[type] = v; }

    void traverse(HMMStateVisitor& visitor);

  private:
    size_t length;
    std::string seq;
    std::string desc;
    const Alphabet& abc;
    std::vector<HMMMatchState> m_state;
    std::vector<HMMDeleteState> d_state;
    std::vector<HMMInsertState> i_state;
    Float1dArray null_emit;
    Float1dArray null_transit[NUM_STATE_TYPE];
    double eff_num;

    void init();
};

class HMMParamEstimator : public HMMStateVisitor{
  public:
    HMMParamEstimator(const TraceVector& traces, const Alphabet& abc, const bool& termi_gap_filter=false);

    void set_seq_weight_estimator(SeqWeightEstimator *estimator)
        { seq_weight_estimator = estimator; }
    void set_eff_seq_num_estimator(EffSeqNumEstimator *estimator)
        { eff_seq_num_estimator = estimator; }
    void set_emit_prob_estimator(EmitProbEstimator *estimator)
        { emit_prob_estimator = estimator; }
    void set_transit_prob_estimator(TransitProbEstimator *estimator)
        { transit_prob_estimator = estimator; }
    void set_msa_filter(MSAFilter *filter)
        { msa_filter = filter; }

    virtual void visit_match(HMMMatchState* state, const size_t& idx);
    virtual void visit_delete(HMMDeleteState* state, const size_t& idx);
    virtual void visit_insert(HMMInsertState* state, const size_t& idx);

    void parameterize(ProfileHMM& model);

  private:
    TraceVector traces;
    const Alphabet& abc;
    const SeqWeightEstimator *seq_weight_estimator;
    const EffSeqNumEstimator *eff_seq_num_estimator;
    const EmitProbEstimator *emit_prob_estimator;
    const TransitProbEstimator *transit_prob_estimator;
    const MSAFilter *msa_filter;
    bool termi_gap_filter;   // flag for terminal gap filtering

    inline Float1dArray calc_seq_weight(const std::vector<std::string>& msa);
    inline double calc_eff_seq_num(const std::vector<std::string>& msa);
    inline Float1dArray calc_emit_freq(const StateType& type, const size_t& idx, const TraceVector& trs, const Float1dArray& wt, const double& eff_num);
    inline Float1dArray calc_tr_freq(const StateType& type, const size_t& idx, const TraceVector& trs, const Float1dArray& wt, const double& eff_num);
    inline Float1dArray estimate_emit_prob(const Float1dArray& f);
    inline Float1dArray estimate_transit_prob(const Float1dArray& f, const StateType& type);
    inline void collect_trace(const StateType type, const size_t& idx, TraceVector& trs, std::vector<std::string>& msa);
};

class PseudocountMixer : public HMMStateVisitor {
  public:
    PseudocountMixer(EmitProbEstimator *emit_prob_estimator,
                     TransitProbEstimator *transit_prob_estimator)
    : emit_prob_estimator(emit_prob_estimator),
      transit_prob_estimator(transit_prob_estimator) {};

    virtual void visit_match(HMMMatchState* state, const size_t& idx);
    virtual void visit_delete(HMMDeleteState* state, const size_t& idx);
    virtual void visit_insert(HMMInsertState* state, const size_t& idx);

    void mix(ProfileHMM& model) { model.traverse(*this); }

  private:
    const EmitProbEstimator *emit_prob_estimator;
    const TransitProbEstimator *transit_prob_estimator;
};

#endif
