#ifndef _ANALYZE_H_
#define _ANALYZE_H_

#include <string>

#include "core.h"
#include "mrf.h"
#include "seq/trace.h"

using std::string;

class MRFModelAnalyzer {
  public:
    MRFModelAnalyzer(const Alphabet& abc) : abc(abc) {}

    //TODO: int calc_pos_score(const string& mrf_filename) { return 0; }
    //TODO: int calc_pair_score(const string& mrf_filename) { return 0; }
    int infer(const string& mrf_filename, const string& seq_filename);

  protected:
    const Alphabet& abc;

    MRF read_mrf(const string& filename);
    TraceVector read_traces(const string& filename);
    double calc_pll(const MRF& model, const string& aseq);
};

#endif
