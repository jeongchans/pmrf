#ifndef _CORE_H_
#define _CORE_H_

#include <fstream>
#include <string>
#include <memory>

#include "seq/alphabet.h"
#include "seq/trace.h"
#include "command.h"
#include "option.h"
#include "mrf.h"
#include "parameterize.h"

using std::string;
using std::shared_ptr;

const AminoAcid AA("-", false, true);
//const AminoAcid AA("=-^", false, true);
//const char GAP_OPEN_SYMBOL = '=';
//const char GAP_EXT_SYMBOL = '-';
//const char GAP_UNALI_SYMBOL = '^';

void assert_file_handle(std::ifstream& ifs, const std::string& filename);
void assert_file_handle(std::ofstream& ofs, const std::string& filename);

class MRFCommandLine;

class MRFCmdProcessor {
  public:
    int run();

  protected:
    MRFCommandLine *cmd_line;
};

class MRFMainProcessor : public MRFCmdProcessor {
  public:
    MRFMainProcessor(int argc, char** argv);
    ~MRFMainProcessor() { delete cmd_line; }

    int run_mrf_cmd(const SubCommand& cmd);

  private:
    int argc;
    char** argv;
};

class MRFBuildProcessor : public MRFCmdProcessor {
  public:
    MRFBuildProcessor(int argc, char** argv);
    ~MRFBuildProcessor() { delete cmd_line; }

    int build();

  private:
    const Alphabet& abc;
    Build::Option* opt;

    shared_ptr<MRFParameterizer> mrf_parameterizer;

    TraceVector read_traces();
    MRF build_mrf(const TraceVector& traces);
    void export_mrf(const MRF& model);
};

class MRFStatProcessor : public MRFCmdProcessor {
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
    MRFStatProcessor(int argc, char** argv);
    ~MRFStatProcessor() { delete cmd_line; }

    int stat();

  private:
    const Alphabet& abc;
    Stat::Option* opt;

    PairScoreVector calc_pair_score(const MRF& model);
    PairScoreVector correct_apc_pair_score(const PairScoreVector& scores);
    PairScoreVector correct_ncps_pair_score(const PairScoreVector& scores);
    void calc_zscore(PairScoreVector& scores);
    vector<FloatType> calc_zscore(const vector<FloatType> scores);
    PosScoreVector calc_pos_score(const MRF& model);
    void calc_zscore(PosScoreVector& scores);
};

class MRFInferProcessor : public MRFCmdProcessor {
  public:
    struct Score {
        double mrf_pll;
        double prof_ll;
    };

    MRFInferProcessor(int argc, char** argv);
    ~MRFInferProcessor() { delete cmd_line; }

    int infer();

  private:
    const Alphabet& abc;
    Infer::Option* opt;

    double calc_mrf_pll(const MRF& model, const string& aseq);
    double calc_prof_ll(const MRF& model, const string& aseq);
    Score calc_score(const MRF& model, const string& aseq);
};

class MRFShowProcessor : public MRFCmdProcessor {
  public:
    MRFShowProcessor(int argc, char** argv);
    ~MRFShowProcessor() { delete cmd_line; }

    int show();

  private:
    const Alphabet& abc;
    Show::Option* opt;
};

#endif
