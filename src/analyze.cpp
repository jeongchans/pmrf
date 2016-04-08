#include "analyze.h"

#include <fstream>
#include <map>
#include <numeric>
#include <iomanip>

#include "mrfio.h"
#include "core.h"
#include "seq/seqio.h"

using std::setw;
using std::setprecision;
using std::fixed;
using std::left;
using std::right;
using std::map;

int MRFModelAnalyzer::infer(const string& mrf_filename, const string& seq_filename) {
    MRF model = read_mrf(mrf_filename);
    double wtscore = calc_pll(model, model.get_seq());
    std::ifstream ifs(seq_filename.c_str());
    assert_file_handle(ifs, seq_filename);
    FastaParser parser(ifs);
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
    while (parser.has_next()) {
        SeqRecord r = parser.next();
        double mtscore = calc_pll(model, r.seq);
        std::cout << setw(4) << right << ++idx << delim
                  << setw(10) << fixed << setprecision(2) << mtscore << delim
                  << setw(10) << fixed << setprecision(2) << mtscore - wtscore << delim
                  << setw(30) << left << r.desc.substr(0, 30) << std::endl;
    }
    return 0;
}

int MRFModelAnalyzer::stat(const string& mrf_filename) {
    string delim = " ";
    MRF model = read_mrf(mrf_filename);
    if (stat_opt.mode == Stat::MODE_PAIR) {
        PairScoreVector scores = calc_pair_score(model);
        if (stat_opt.corr == Stat::CORR_APC) scores = correct_apc_pair_score(scores);
        else if (stat_opt.corr == Stat::CORR_NCPS) scores = correct_ncps_pair_score(scores);
        calc_zscore(scores);
        std::cout << setw(4) << right << "Pos1" << delim
                  << setw(4) << "Pos2" << delim
                  << setw(10) << "Score" << delim
                  << setw(10) << "Z-score" << std::endl;
        std::cout << "----" << delim
                  << "----" << delim
                  << "----------" << delim
                  << "----------" << std::endl;
        for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
            std::cout << setw(4) << right << pos->idx.idx1 + 1 << delim
                      << setw(4) << pos->idx.idx2 + 1 << delim
                      << setw(10) << fixed << pos->score << delim
                      << setw(10) << fixed << pos->zscore << std::endl;
        }
    } else if (stat_opt.mode == Stat::MODE_POS) {
        PosScoreVector scores = calc_pos_score(model);
        calc_zscore(scores);
        std::cout << setw(4) << right << "Pos" << delim
                  << setw(10) << "Score" << delim
                  << setw(10) << "Z-score" << std::endl;
        std::cout << "----" << delim
                  << "----------" << delim
                  << "----------" << std::endl;
        for (PosScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
            std::cout << setw(4) << right << pos->idx + 1 << delim
                      << setw(10) << fixed << pos->score << delim
                      << setw(10) << fixed << pos->zscore << std::endl;
        }
    }
    return 0;
}

MRF MRFModelAnalyzer::read_mrf(const string& filename) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return MRFImporter().import_model(ifs, abc);
}

double MRFModelAnalyzer::calc_pll(const MRF& model, const string& aseq) {
    double pll = 0.;
    size_t n = model.get_length();
    const Alphabet& abc = model.get_alphabet();
    for (size_t i = 0; i < n; ++i) {
        FloatType w1;
        string s1 = abc.get_degeneracy(aseq[i], &w1);
        for (string::iterator c = s1.begin(); c != s1.end(); ++c)
            pll += w1 * model.get_node(i).get_weight()(abc.get_idx(*c));
        for (size_t j = i + 1; j < n; ++j) {
            if (model.has_edge(i, j)) {
                FloatType w2;
                string s2 = abc.get_degeneracy(aseq[j], &w2);
                for (string::iterator c1 = s1.begin(); c1 != s1.end(); ++c1) {
                    for (string::iterator c2 = s2.begin(); c2 != s2.end(); ++c2)
                        pll += w1 * w2 * model.get_edge(i, j).get_weight()(abc.get_idx(*c1), abc.get_idx(*c2));
                }
            }
        }
    }
    return pll;
}

MRFModelAnalyzer::PairScoreVector MRFModelAnalyzer::calc_pair_score(const MRF& model) {
    PairScoreVector scores;
    size_t n = model.get_length();
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (model.has_edge(i, j)) {
                const MatrixXf& w = model.get_edge(i, j).get_weight();
                scores.push_back(PairScore(EdgeIndex(i, j), w.topLeftCorner<20, 20>().norm()));
            }
        }
    }
    return scores;
}

MRFModelAnalyzer::PairScoreVector MRFModelAnalyzer::correct_apc_pair_score(const PairScoreVector& scores) {
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

MRFModelAnalyzer::PairScoreVector MRFModelAnalyzer::correct_ncps_pair_score(const PairScoreVector& scores) {
    map<size_t, map<size_t, FloatType> > cm;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
        cm[pos->idx.idx1][pos->idx.idx2] = pos->score;
        cm[pos->idx.idx2][pos->idx.idx1] = pos->score;
    }
    PairScoreVector ncps;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) {
        map<size_t, FloatType>& cm1 = cm[pos->idx.idx1];
        map<size_t, FloatType>& cm2 = cm[pos->idx.idx2];
        FloatType s = 0.;
        size_t n = 0;
        for (map<size_t, FloatType>::const_iterator p = cm1.begin(); p != cm1.end(); ++p) {
            if (p->first != pos->idx.idx2 && cm2.find(p->first) != cm2.end()) {
                s += p->second * cm2[p->first];
                ++n;
            }
        }
        if (n > 0) s /= n;
        ncps.push_back(PairScore(pos->idx, s));
    }
    FloatType tot_mn = 0.;
    for (PairScoreVector::const_iterator pos = ncps.begin(); pos != ncps.end(); ++pos) tot_mn += pos->score;
    tot_mn /= ncps.size();
    tot_mn = sqrt(tot_mn);
    for (PairScoreVector::iterator pos = ncps.begin(); pos != ncps.end(); ++pos) pos->score /= tot_mn;
    PairScoreVector ret;
    for (size_t i = 0; i < scores.size(); ++i)
        ret.push_back(PairScore(scores[i].idx, scores[i].score - ncps[i].score));
    return ret;
}

void MRFModelAnalyzer::calc_zscore(PairScoreVector& scores) {
    vector<FloatType> vec;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) vec.push_back(pos->score);
    vec = calc_zscore(vec);
    for (size_t i = 0; i < scores.size(); ++i) scores[i].zscore = vec[i];
}

vector<FloatType> MRFModelAnalyzer::calc_zscore(const vector<FloatType> scores) {
    FloatType mn = accumulate(scores.begin(), scores.end(), 0.) / scores.size();
    FloatType mn2 = inner_product(scores.begin(), scores.end(), scores.begin(), 0.) / scores.size();
    FloatType sd = sqrt(mn2 - mn * mn);
    vector<FloatType> ret;
    for (vector<FloatType>::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) ret.push_back((*pos - mn) / sd);
    return ret;
}

MRFModelAnalyzer::PosScoreVector MRFModelAnalyzer::calc_pos_score(const MRF& model) {
    PosScoreVector scores;
    size_t n = model.get_length();
    const MatrixXf psfm = model.get_psfm();
    for (size_t i = 0; i < n; ++i) {
        const VectorXf& w = model.get_node(i).get_weight().head<20>();
        VectorXf p = w.unaryExpr(&exp);
        p /= p.sum();
        VectorXf f = psfm.row(i).head<20>();
        scores.push_back(PosScore(i, kld(f, p)));
    }
    return scores;
}

void MRFModelAnalyzer::calc_zscore(PosScoreVector& scores) {
    vector<FloatType> vec;
    for (PosScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) vec.push_back(pos->score);
    vec = calc_zscore(vec);
    for (size_t i = 0; i < scores.size(); ++i) scores[i].zscore = vec[i];
}
