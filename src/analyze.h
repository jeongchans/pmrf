#ifndef _ANALYZE_H_
#define _ANALYZE_H_

#include <string>

#include "core.h"
#include "mrf.h"
#include "seq/trace.h"

using std::string;

enum StatMode { STATMODE_PAIR, STATMODE_POS };
enum StatCorrect { STATCORR_NONE, STATCORR_APC, STATCORR_NCPS };

class MRFModelAnalyzer {

  public:
    class StatOption {
      public:
        StatOption(const StatMode& mode=STATMODE_PAIR, const StatCorrect& corr=STATCORR_APC, const bool& zscore=true)
        : mode(mode), corr(corr), zscore(zscore) {};
        
        StatMode mode;
        StatCorrect corr;
        bool zscore;
    };

    struct PairScore {
        EdgeIndex idx;
        FloatType score;

        PairScore(const EdgeIndex& idx, const FloatType& score) : idx(idx), score(score) {}
    };
    typedef vector<PairScore> PairScoreVector;

    struct PosScore {
        size_t idx;
        FloatType score;

        PosScore(const size_t& idx, const FloatType& score) : idx(idx), score(score) {}
    };
    typedef vector<PosScore> PosScoreVector;

  public:
    MRFModelAnalyzer(const Alphabet& abc) : abc(abc) {}

    int infer(const string& mrf_filename, const string& seq_filename);
    int stat(const string& mrf_filename);

    StatOption stat_opt;

  protected:
    const Alphabet& abc;

    MRF read_mrf(const string& filename);
    TraceVector read_traces(const string& filename);
    double calc_pll(const MRF& model, const string& aseq);
    PairScoreVector calc_pair_score(const MRF& model);
    PairScoreVector correct_pair_score(const PairScoreVector& scores);
    vector<FloatType> calc_zscore(const vector<FloatType> scores);
    PosScoreVector calc_pos_score(const MRF& model);
};

#endif
