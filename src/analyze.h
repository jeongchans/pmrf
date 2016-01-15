#ifndef _ANALYZE_H_
#define _ANALYZE_H_

#include <string>

#include "core.h"

using std::string;

class MRFModelAnalyzer {
  public:
    MRFModelAnalyzer(const Alphabet& abc) : abc(abc) {}

    //TODO: int calc_pos_score(const string& mrf_filename) { return 0; }
    //TODO: int calc_pair_score(const string& mrf_filename) { return 0; }
    int infer(const string& mrf_filename, const string& seq_filename);

  private:
    const Alphabet& abc;
};

#endif
