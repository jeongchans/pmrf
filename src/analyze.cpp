#include "analyze.h"

#include <fstream>
#include <map>
#include <numeric>

#include "mrfio.h"
#include "core.h"

using std::setw;
using std::setprecision;
using std::fixed;
using std::left;
using std::right;
using std::map;

int MRFModelAnalyzer::infer(const string& mrf_filename, const string& seq_filename) {
    MRF model = read_mrf(mrf_filename);
    double wtscore = calc_pll(model, model.get_seq());
    TraceVector traces = read_traces(seq_filename);
    string delim = " ";
    std::cout << setw(4) << right << "#" << delim
              << setw(10) << "Score" << delim
              << setw(10) << "Diff" << delim
              << setw(30) << left << "Description" << std::endl;
    std::cout << "----" << delim
              << "----------" << delim
              << "----------" << delim
              << "--------------------" << std::endl;
    size_t idx = 0;
    for (TraceVector::const_iterator pos = traces.begin(); pos != traces.end(); ++pos) {
        double mtscore = calc_pll(model, pos->get_matched_aseq());
        std::cout << setw(4) << right << ++idx << delim
                  << setw(10) << fixed << setprecision(2) << mtscore << delim
                  << setw(10) << fixed << setprecision(2) << mtscore - wtscore << delim
                  << setw(30) << left << pos->get_desc().substr(0, 30) << std::endl;
    }
    return 0;
}

int MRFModelAnalyzer::stat_pair(const string& mrf_filename) {
    MRF model = read_mrf(mrf_filename);
    PairScoreVector scores = calc_pair_score(model);
    scores = correct_pair_score(scores);
    vector<FloatType> vec;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) vec.push_back(pos->score);
    vec = calc_zscore(vec);
    for (size_t i = 0; i < scores.size(); ++i) scores[i].score = vec[i];
    string delim = " ";
    std::cout << setw(4) << right << "Idx1" << delim
              << setw(4) << "Idx2" << delim
              << setw(10) << "Score" << std::endl;
    std::cout << "----" << delim
              << "----" << delim
              << "----------" << std::endl;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
        std::cout << setw(4) << right << pos->idx.idx1 + 1 << delim
                  << setw(4) << pos->idx.idx2 + 1 << delim
                  << setw(10) << fixed << pos->score << std::endl;
    }
    return 0;
}

MRF MRFModelAnalyzer::read_mrf(const string& filename) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return MRFImporter().import_model(ifs, abc);
}

TraceVector MRFModelAnalyzer::read_traces(const string& filename) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return TraceImporter(abc).import(ifs, AFASTA).second;
}

double MRFModelAnalyzer::calc_pll(const MRF& model, const string& aseq) {
    double pll = 0.;
    size_t n = model.get_length();
    for (size_t i = 0; i < n; ++i) {
        int idx1 = model.get_var_idx(aseq[i]);
        pll += model.get_node(i).get_weight()(idx1);
        for (size_t j = i + 1; j < n; ++j) {
            int idx2 = model.get_var_idx(aseq[j]);
            if (model.has_edge(i, j)) pll += model.get_edge(i, j).get_weight()(idx1, idx2);
        }
    }
    return pll;
}

MRFModelAnalyzer::PairScoreVector MRFModelAnalyzer::calc_pair_score(const MRF& model) {
    PairScoreVector scores;
    size_t n = model.get_length();
    Range r(blitz::fromStart, 19);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (model.has_edge(i, j)) {
                const Float2dArray& w = model.get_edge(i, j).get_weight();
                scores.push_back(PairScore(EdgeIndex(i, j), norm(w(r, r))));
            }
        }
    }
    return scores;
}

MRFModelAnalyzer::PairScoreVector MRFModelAnalyzer::correct_pair_score(const PairScoreVector& scores) {
    map<size_t, FloatType> mn;
    map<size_t, size_t> num_edge;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
        mn[pos->idx.idx1] = 0.;
        num_edge[pos->idx.idx1] = 0;
        mn[pos->idx.idx2] = 0.;
        num_edge[pos->idx.idx2] = 0;
    }
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
        mn[pos->idx.idx1] += pos->score;
        num_edge[pos->idx.idx1]++;
        mn[pos->idx.idx2] += pos->score;
        num_edge[pos->idx.idx2]++;
    }
    for (map<size_t, FloatType>::iterator pos = mn.begin(); pos != mn.end(); ++pos) pos->second /= num_edge[pos->first];
    FloatType tot_mn = 0.;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) tot_mn += pos->score;
    tot_mn /= scores.size();
    PairScoreVector ret;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
        ret.push_back(PairScore(pos->idx, pos->score - mn[pos->idx.idx1] * mn[pos->idx.idx2] / tot_mn));
    }
    return ret;
}

vector<FloatType> MRFModelAnalyzer::calc_zscore(const vector<FloatType> scores) {
    FloatType mn = accumulate(scores.begin(), scores.end(), 0.) / scores.size();
    FloatType mn2 = inner_product(scores.begin(), scores.end(), scores.begin(), 0.) / scores.size();
    FloatType sd = sqrt(mn2 - mn * mn);
    vector<FloatType> ret;
    for (vector<FloatType>::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) ret.push_back((*pos - mn) / sd);
    return ret;
}
