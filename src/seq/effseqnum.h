#ifndef _EFFSEQNUM_H_
#define _EFFSEQNUM_H_

#include "seq/alphabet.h"
#include "seqweight.h"

class EffSeqNumEstimator {
  public:
    virtual double estimate(const vector<string>& msa) const = 0;
    virtual double estimate(const vector<string>& msa, const VectorXf& wt) const = 0;
};

// Null effective sequence number estimator gives the total number of MSA sequences
class NullEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    NullEffSeqNumEstimator() {};

    virtual double estimate(const vector<string>& msa) const { return msa.size(); }
    virtual double estimate(const vector<string>& msa, const VectorXf& wt) const { return estimate(msa); }
};

// Average number of different residue types is used as effective sequence
// number measure
class RTEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    RTEffSeqNumEstimator(const Alphabet& abc) : abc(abc) {};

    virtual double estimate(const vector<string>& msa) const;
    virtual double estimate(const vector<string>& msa, const VectorXf& wt) const { return estimate(msa); }

  private:
    const Alphabet& abc;

    size_t calc_num_res_type(const string& column) const;
};

// Exponential of average entropy
// Biegert and Soding. Proc Natl Acad Sci USA (2009) vol. 106 (10) pp. 3770-5
// Note: By changing the base of logarithm as e, the resulting Neff is bounded
//       to [1, 20].
//       That is, the Neff indicates the expected number of different amino acids
//       observed at a position.
class ExpEntropyEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    ExpEntropyEffSeqNumEstimator(const Alphabet& abc, SeqWeightEstimator* seq_weight_estimator) : abc(abc), seq_weight_estimator(seq_weight_estimator) {};

    virtual double estimate(const vector<string>& msa) const { return estimate(msa, seq_weight_estimator->estimate(msa)); }
    virtual double estimate(const vector<string>& msa, const VectorXf& wt) const;

  private:
    const Alphabet& abc;
    SeqWeightEstimator *seq_weight_estimator;
};

// Exponential of average joint entropy
// Note: The resulting Neff is bounded to [1, 400].
//       That is, the Neff indicates the expected number of different amino acid pairs
//       observed at a position pair.
class ExpJointEntropyEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    ExpJointEntropyEffSeqNumEstimator(const Alphabet& abc, SeqWeightEstimator* seq_weight_estimator) : abc(abc), seq_weight_estimator(seq_weight_estimator) {};

    virtual double estimate(const vector<string>& msa) const { return estimate(msa, seq_weight_estimator->estimate(msa)); }
    virtual double estimate(const vector<string>& msa, const VectorXf& wt) const;

  private:
    const Alphabet& abc;
    SeqWeightEstimator *seq_weight_estimator;
};

// Number of sequences in the clustered MSA by a given sequence identity
// Note: The maxidt specifies the maximum sequence identity in the clustered MSA
class ClstrEffSeqNumEstimator : public EffSeqNumEstimator {
  public:
    ClstrEffSeqNumEstimator(const Alphabet& abc, const float maxidt) : abc(abc), maxidt(maxidt) {};

    virtual double estimate(const vector<string>& msa) const;
    virtual double estimate(const vector<string>& msa, const VectorXf& wt) const { return estimate(msa); }

  private:
    const Alphabet& abc;
    float maxidt;       // maximum sequence identity
};

#endif
