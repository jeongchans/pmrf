#ifndef _ANALYZE_H_
#define _ANALYZE_H_

#include <string>

#include "core.h"
#include "mrf.h"
#include "seq/trace.h"

using std::string;

class MRFModelAnalyzer {

  public:
    struct PairScore {
        EdgeIndex idx;
        FloatType score;

        PairScore(const EdgeIndex& idx, const FloatType& score) : idx(idx), score(score) {}
    };
    typedef vector<PairScore> PairScoreVector;

  public:
    MRFModelAnalyzer(const Alphabet& abc) : abc(abc) {}

    int infer(const string& mrf_filename, const string& seq_filename);
    int stat_pair(const string& mrf_filename);
    //TODO: int calc_pos_score(const string& mrf_filename) { return 0; }

  protected:
    const Alphabet& abc;

    MRF read_mrf(const string& filename);
    TraceVector read_traces(const string& filename);
    double calc_pll(const MRF& model, const string& aseq);
    PairScoreVector calc_pair_score(const MRF& model);
    PairScoreVector correct_pair_score(const PairScoreVector& scores);
    vector<FloatType> calc_zscore(const vector<FloatType> scores);
};

#endif
