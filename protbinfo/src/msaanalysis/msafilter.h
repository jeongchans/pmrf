#ifndef _MSAFILTER_H_
#define _MSAFILTER_H_

#include "common/common.h"
#include "seq/alphabet.h"

class MSAFilter {
  public:
    virtual vector<string> filter(const vector<string>& msa) const = 0;
};

// Null MSA filter does nothing
class NullMSAFilter : public MSAFilter {
  public:
    NullMSAFilter() {};

    virtual vector<string> filter(const vector<string>& msa) const { return msa; }
};

// Terminal gap remover
// This MSA filter removes terminal residues by maximum terminal gap percentage
class TerminalGapRemover : public MSAFilter {
  public:
    TerminalGapRemover(const Alphabet& abc, const double& max_gap_perc) 
    : abc(abc), max_gap_perc(max_gap_perc), TERMI('T'), NON_TERMI('S') {};

    virtual vector<string> filter(const vector<string>& msa) const;

  private:
    const Alphabet& abc;
    const double max_gap_perc;
    const char TERMI;         // terminal gap
    const char NON_TERMI;     // match or delete

    vector<string> build_terminus_state(const vector<string>& msa) const;
};

#endif
