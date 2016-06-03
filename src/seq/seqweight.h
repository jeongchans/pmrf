#ifndef _SEQWEIGHT_H_
#define _SEQWEIGHT_H_

#include "util/common.h"
#include "util/numeric.h"
#include "seq/alphabet.h"

class SeqWeightEstimator {
  public:
    virtual VectorXf estimate(const vector<string>& msa) const = 0;
};

// Null sequence weight estimator gives the same weights for all MSA sequences
class NullSeqWeightEstimator : public SeqWeightEstimator {
  public:
    NullSeqWeightEstimator() {};

    virtual VectorXf estimate(const vector<string>& msa) const {
        size_t n = msa.size();
        return VectorXf::Constant(n, 1. / n);
    }
};

// Henikoff and Henikoff's position based sequence weights
// JMB-243-574 (1994)
//
// The original algorithm does not consider gap and non-canonical symbols.
// To handle them, we ignore those symbols, refering HMMER implementation.
// This treatment can increase the weights for logner sequences. Thus,
// sequence weights are additionally normalized with the sequence length.
class PBSeqWeightEstimator : public SeqWeightEstimator {
  public:
    PBSeqWeightEstimator() : dim(26) {};

    virtual VectorXf estimate(const vector<string>& msa) const;

  private:
    const int dim;

    VectorXf calc_residue_weight(const vector<string>& msa, const int& idx) const;
    bool is_allowed(const char& c) const;
    int abc_idx(const char& c) const;
};

// Modified version of Henikoff and Henikoff's position based sequence weights
// NAR-25-3389 (1997)
//
// PSI-BLAST's sequence weight scheme. This algorithm treats gap as an
// additional symbol, and ignores any column consisting of an identical 
// residue.
class PSIBLASTPBSeqWeightEstimator : public SeqWeightEstimator {
  public:
    PSIBLASTPBSeqWeightEstimator() : dim(27), gap_idx(26) {};

    virtual VectorXf estimate(const vector<string>& msa) const;

  private:
    const int dim;
    const int gap_idx;

    VectorXf calc_residue_weight(const vector<string>& msa, const int& idx) const;
    int abc_idx(const char& c) const;
};

// Number of sequences in the clustered MSA by a given sequence identity
// Note: The maxidt specifies the maximum sequence identity in the clustered MSA
class ClstrSeqWeightEstimator : public SeqWeightEstimator {
  public:
    ClstrSeqWeightEstimator(const Alphabet& abc, const float maxidt) : abc(abc), maxidt(maxidt) {};

    virtual VectorXf estimate(const vector<string>& msa) const;

  private:
    const Alphabet& abc;
    float maxidt;       // maximum sequence identity

    float calc_identity(const string& seq1, const string& seq2) const;
};

#endif
