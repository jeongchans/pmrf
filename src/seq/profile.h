#ifndef _SEQ_PROFILE_H_
#define _SEQ_PROFILE_H_

#include "util/common.h"
#include "trace.h"
#include "alphabet.h"
#include "seqweight.h"
#include "effseqnum.h"
#include "msafilter.h"
#include "probdistrib.h"

class Profile {
  public:
    Profile(const size_t& length, const Alphabet& abc);
    Profile(const string& seq, const Alphabet& abc);

    const Float2dArray get_prob() const { return prob; }

    const Float1dArray get_prob(const size_t& idx) const { return prob(idx, ALL); }
    void set_prob(const size_t& idx, const Float1dArray& v) { prob(idx, ALL) = v; }

    double get_eff_num(const size_t& idx) const { return eff_num(idx); }
    void set_eff_num(const size_t& idx, const double& neff) { eff_num(idx) = neff; }

    size_t get_length() const { return length; }
    string get_seq() const { return seq; }

  private:
    size_t length;
    string seq;
    const Alphabet& abc;
    Float2dArray prob;
    Float1dArray eff_num;

    void init();
};

class ProfileBuilder {
  public:
    ProfileBuilder(const Alphabet& abc);

    Profile build(const TraceVector& traces) const;

    void set_seq_weight_estimator(SeqWeightEstimator *estimator)
        { seq_weight_estimator = estimator; }
    void set_eff_seq_num_estimator(EffSeqNumEstimator *estimator)
        { eff_seq_num_estimator = estimator; }
    void set_msa_filter(MSAFilter *filter)
        { msa_filter = filter; }
    void set_emit_prob_estimator(EmitProbEstimator *estimator)
        { emit_prob_estimator = estimator; }

  private:
    const Alphabet& abc;
    const SeqWeightEstimator *seq_weight_estimator;
    const EffSeqNumEstimator *eff_seq_num_estimator;
    const MSAFilter *msa_filter;
    const EmitProbEstimator *emit_prob_estimator;

    void collect_trace(const TraceVector& traces, const size_t& idx, TraceVector& r_trs, vector<string>& r_msa) const;
    Float1dArray calc_freq(const TraceVector& trs, const size_t& idx, const Float1dArray& wt, const double& eff_num) const;
};

#endif
