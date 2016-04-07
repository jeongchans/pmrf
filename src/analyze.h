#ifndef _ANALYZE_H_
#define _ANALYZE_H_

#include <string>

#include "core.h"
#include "option.h"
#include "mrf.h"
#include "seq/trace.h"

using std::string;

class MRFModelAnalyzer {

  public:
    struct PairScore {
        EdgeIndex idx;
        FloatType score;
        FloatType zscore;

        PairScore(const EdgeIndex& idx, const FloatType& score) : idx(idx), score(score) {}
    };
    typedef vector<PairScore> PairScoreVector;

    struct PosScore {
        size_t idx;
        FloatType score;
        FloatType zscore;

        PosScore(const size_t& idx, const FloatType& score) : idx(idx), score(score) {}
    };
    typedef vector<PosScore> PosScoreVector;

  public:
    MRFModelAnalyzer(const Alphabet& abc) : abc(abc) {}

    int infer(const string& mrf_filename, const string& seq_filename);
    int stat(const string& mrf_filename);

    Stat::Option stat_opt;

  protected:
    const Alphabet& abc;

    MRF read_mrf(const string& filename);
    double calc_pll(const MRF& model, const string& aseq);
    PairScoreVector calc_pair_score(const MRF& model);
    PairScoreVector correct_apc_pair_score(const PairScoreVector& scores);
    PairScoreVector correct_ncps_pair_score(const PairScoreVector& scores);
    void calc_zscore(PairScoreVector& scores);
    vector<FloatType> calc_zscore(const vector<FloatType> scores);
    PosScoreVector calc_pos_score(const MRF& model);
    void calc_zscore(PosScoreVector& scores);
};

#endif
