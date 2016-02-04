#ifndef _EFFSEQNUM_H_
#define _EFFSEQNUM_H_

#include "seq/alphabet.h"
#include "seqweight.h"

class EffSeqNumEstimator {
  public:
    virtual double estimate(const vector<string>& msa) const = 0;
};

// Null effective sequence number estimator gives the total number of MSA sequences
class NullEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    NullEffSeqNumEstimator() {};

    virtual double estimate(const vector<string>& msa) const {
        return msa.size();
    }
};

// Average number of different residue types is used as effective sequence
// number measure
class RTEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    RTEffSeqNumEstimator(const Alphabet& abc) : abc(abc) {};

    virtual double estimate(const vector<string>& msa) const;

  private:
    const Alphabet& abc;

    size_t calc_num_res_type(const string& column) const;
};

// Exponential of average entropy
// Biegert and Soding. Proc Natl Acad Sci USA (2009) vol. 106 (10) pp. 3770-5
// Note: By changing the base of logarithm as e, the resulting Neff is bounded
//       to [1, 20].
class ExpEntropyEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    ExpEntropyEffSeqNumEstimator(const Alphabet& abc, SeqWeightEstimator* seq_weight_estimator) : abc(abc), seq_weight_estimator(seq_weight_estimator) {};

    virtual double estimate(const vector<string>& msa) const;

  private:
    const Alphabet& abc;
    SeqWeightEstimator *seq_weight_estimator;
};

#endif