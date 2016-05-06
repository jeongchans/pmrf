#include "core.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <map>
#include <numeric>
#include <iomanip>

#include "mrfio.h"
#include "seq/seqio.h"

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::setprecision;
using std::fixed;
using std::left;
using std::right;
using std::map;

MRF read_mrf(const string& filename, const Alphabet& abc) {
    std::ifstream ifs(filename.c_str());
    assert_file_handle(ifs, filename);
    return MRFImporter().import_model(ifs, abc);
}

int MRFCmdProcessor::run() {
    if (cmd_line->is_valid()) {
        return cmd_line->process_command(this);
    } else {
        cmd_line->show_error();
        exit(EXIT_FAILURE);
    }
}

void assert_file_handle(std::ifstream& ifs, const std::string& filename) {
    if (!ifs) {
        cerr << "can't open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
}

void assert_file_handle(std::ofstream& ofs, const std::string& filename) {
    if (!ofs) {
        cerr << "can't open file: " << filename << endl;
        exit(EXIT_FAILURE);
    }
}

/**
   @class MRFMainProcessor
 */

MRFMainProcessor::MRFMainProcessor(int argc, char** argv) {
    cmd_line = new MRFMainCommandLine(argc, argv);
    this->argc = argc;
    this->argv = argv;
}

int MRFMainProcessor::run_mrf_cmd(const SubCommand& cmd) {
    argv++;
    argc--;
    switch (cmd) {
    case NONE:
    case HELP:
        cmd_line->show_help();
        exit(0);
    case BUILD:
        return MRFBuildProcessor(argc, argv).run();
    case INFER:
        return MRFInferProcessor(argc, argv).run();
    case STAT:
        return MRFStatProcessor(argc, argv).run();
    case SHOW:
        return MRFShowProcessor(argc, argv).run();
    }
    return 1;
}

/**
   @class MRFBuildProcessor
 */

MRFBuildProcessor::MRFBuildProcessor(int argc, char** argv) : abc(AA) {
    cmd_line = new MRFBuildCommandLine(argc, argv);
    opt = &(((MRFBuildCommandLine*)cmd_line)->opt);
}

int MRFBuildProcessor::build() {
    TraceVector traces = read_traces();
    MRF model = build_mrf(traces);
    export_mrf(model);
    return 0;
}

TraceVector MRFBuildProcessor::read_traces() {
    std::ifstream ifs(opt->msa_filename.c_str());
    assert_file_handle(ifs, opt->msa_filename);
    return TraceImporter(abc).import(ifs, opt->msa_fmt).second;
}

MRF MRFBuildProcessor::build_mrf(const TraceVector& traces) {
    EdgeIndexVector eidxs;
    EdgeIndexVector* p=NULL;
    if (!opt->eidx_filename.empty()) {
        std::ifstream ifs(opt->eidx_filename.c_str());
        assert_file_handle(ifs, opt->eidx_filename);
        eidxs = EdgeIndexImporter().import(ifs);
        p = &eidxs;
    }
    MSAAnalyzer msa_analyzer(opt->msa_analyzer_opt, abc);
    MRF model(traces[0].get_seq(), abc, p);
    MRFParameterizer mrf_parameterizer(msa_analyzer);
    mrf_parameterizer.opt = opt->parameterizer_opt;
    mrf_parameterizer.optim_opt = opt->optim_opt;
    mrf_parameterizer.parameterize(model, traces);
    return model;
}

void MRFBuildProcessor::export_mrf(const MRF& model) {
    if (opt->out_filename.empty()) cout << model;
    else {
        MRFExporter mrf_exporter;
        std::ofstream ofs(opt->out_filename.c_str());
        assert_file_handle(ofs, opt->out_filename);
        mrf_exporter.export_model(model, ofs);
    }
}

/**
   @class MRFStatProcessor
 */

MRFStatProcessor::MRFStatProcessor(int argc, char** argv) : abc(AA) {
    cmd_line = new MRFStatCommandLine(argc, argv);
    opt = &(((MRFStatCommandLine*)cmd_line)->opt);
}

int MRFStatProcessor::stat() {
    string delim = " ";
    MRF model = read_mrf(opt->mrf_filename, abc);
    if (opt->mode == Stat::MODE_PAIR) {
        PairScoreVector scores = calc_pair_score(model);
        if (opt->corr == Stat::CORR_APC) scores = correct_apc_pair_score(scores);
        else if (opt->corr == Stat::CORR_NCPS) scores = correct_ncps_pair_score(scores);
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
    } else if (opt->mode == Stat::MODE_POS) {
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

MRFStatProcessor::PairScoreVector MRFStatProcessor::calc_pair_score(const MRF& model) {
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

MRFStatProcessor::PairScoreVector MRFStatProcessor::correct_apc_pair_score(const PairScoreVector& scores) {
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

MRFStatProcessor::PairScoreVector MRFStatProcessor::correct_ncps_pair_score(const PairScoreVector& scores) {
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

void MRFStatProcessor::calc_zscore(PairScoreVector& scores) {
    vector<FloatType> vec;
    for (PairScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) vec.push_back(pos->score);
    vec = calc_zscore(vec);
    for (size_t i = 0; i < scores.size(); ++i) scores[i].zscore = vec[i];
}

vector<FloatType> MRFStatProcessor::calc_zscore(const vector<FloatType> scores) {
    FloatType mn = accumulate(scores.begin(), scores.end(), 0.) / scores.size();
    FloatType mn2 = inner_product(scores.begin(), scores.end(), scores.begin(), 0.) / scores.size();
    FloatType sd = sqrt(mn2 - mn * mn);
    vector<FloatType> ret;
    for (vector<FloatType>::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) ret.push_back((*pos - mn) / sd);
    return ret;
}

MRFStatProcessor::PosScoreVector MRFStatProcessor::calc_pos_score(const MRF& model) {
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

void MRFStatProcessor::calc_zscore(PosScoreVector& scores) {
    vector<FloatType> vec;
    for (PosScoreVector::const_iterator pos = scores.begin(); pos != scores.end(); ++pos) vec.push_back(pos->score);
    vec = calc_zscore(vec);
    for (size_t i = 0; i < scores.size(); ++i) scores[i].zscore = vec[i];
}

/**
   @class MRFInferProcessor
 */

MRFInferProcessor::MRFInferProcessor(int argc, char** argv) : abc(AA) {
    cmd_line = new MRFInferCommandLine(argc, argv);
    opt = &(((MRFInferCommandLine*)cmd_line)->opt);
}

int MRFInferProcessor::infer() {
    MRF model = read_mrf(opt->mrf_filename, abc);
    Score ref_score = calc_score(model, model.get_seq());
    const string& seq_filename = opt->seq_filename;
    std::ifstream ifs(seq_filename.c_str());
    assert_file_handle(ifs, seq_filename);
    FastaParser parser(ifs);
    string delim = " ";
    static const string header =
        "     -------- MRF -------- ------ Profile ------                     \n"
        "   #      Score       Diff      Score       Diff Description         \n"
        "---- ---------- ---------- ---------- ---------- --------------------\n";
    std::cout << header;
    size_t idx = 0;
    while (parser.has_next()) {
        SeqRecord r = parser.next();
        Score score = calc_score(model, r.seq);
        std::cout << setw(4) << right << ++idx << delim
                  << setw(10) << fixed << setprecision(2) << score.mrf_pll << delim
                  << setw(10) << fixed << setprecision(2) << score.mrf_pll - ref_score.mrf_pll << delim
                  << setw(10) << fixed << setprecision(2) << score.prof_ll << delim
                  << setw(10) << fixed << setprecision(2) << score.prof_ll - ref_score.prof_ll << delim
                  << setw(30) << left << r.desc.substr(0, 30) << std::endl;
    }
    return 0;
}

MRFInferProcessor::Score MRFInferProcessor::calc_score(const MRF& model, const string& aseq) {
    Score score;
    score.mrf_pll = calc_mrf_pll(model, aseq);
    score.prof_ll = calc_prof_ll(model, aseq);
    return score;
}

double MRFInferProcessor::calc_mrf_pll(const MRF& model, const string& aseq) {
    double pll = 0.;
    size_t n = model.get_length();
    const Alphabet& abc = model.get_alphabet();
    for (size_t i = 0; i < n; ++i) {
        FloatType w1;
        const string s1 = abc.get_degeneracy(aseq[i], &w1);
        for (auto c = s1.cbegin(); c != s1.cend(); ++c)
            pll += w1 * (model.get_node(i).get_weight().array() - opt->node_offset)(abc.get_idx(*c));
        for (size_t j = i + 1; j < n; ++j) {
            if (model.has_edge(i, j)) {
                FloatType w2;
                const string s2 = abc.get_degeneracy(aseq[j], &w2);
                for (auto c1 = s1.begin(); c1 != s1.end(); ++c1) {
                    for (auto c2 = s2.begin(); c2 != s2.end(); ++c2)
                        pll += w1 * w2 * (model.get_edge(i, j).get_weight().array() - opt->edge_offset)(abc.get_idx(*c1), abc.get_idx(*c2));
                }
            }
        }
    }
    return pll;
}

double MRFInferProcessor::calc_prof_ll(const MRF& model, const string& aseq) {
    double ll = 0.;
    size_t n = model.get_length();
    const Alphabet& abc = model.get_alphabet();
    MatrixXf mat = model.get_psfm().unaryExpr(&log);
    mat = (mat.array() - opt->prof_offset).matrix();
    mat.col(20) = VectorXf::Constant(mat.rows(), opt->gap_score);
    for (size_t i = 0; i < n; ++i) {
        float w;
        string s = abc.get_degeneracy(aseq[i], &w);
        for (auto c = s.cbegin(); c != s.cend(); ++c) ll += w * mat(i, abc.get_idx(*c));
    }
    return ll;
}

/**
   @class MRFShowProcessor
 */

MRFShowProcessor::MRFShowProcessor(int argc, char** argv) : abc(AA), sep("\t"), prec(4) {
    cmd_line = new MRFShowCommandLine(argc, argv);
    opt = &(((MRFShowCommandLine*)cmd_line)->opt);
}

int MRFShowProcessor::show() {
    MRF model = read_mrf(opt->mrf_filename, abc);
    if (opt->seq_flag) show_seq(cout, model);
    if (opt->prof_flag) show_prof(cout, model);
    if (opt->mrf_flag) show_mrf(cout, model);
    if (opt->node_flag) show_mrf_node(cout, model, opt->v_pos);
    if (opt->edge_flag) show_mrf_edge(cout, model, opt->w_pos);
    if (!opt->seq_flag && !opt->prof_flag && !opt->mrf_flag && !opt->node_flag && !opt->edge_flag) cout << model;
    return 0;
}

void MRFShowProcessor::show_seq(ostream& os, const MRF& model) {
    size_t n = model.get_length();
    string seq = model.get_seq();
    string fseq;
    const size_t block_width = 10;
    const size_t line_width = 50;
    for (size_t i = 0; i < n; i += block_width) {
        fseq += seq.substr(i, block_width) + " ";
        if ((i + block_width) % line_width == 0) {
            string tick = std::to_string(i + block_width);
            while (tick.size() < 4) tick = " " + tick;
            fseq += tick + "\n";
        }
    }
    os << "# SEQUENCE" << endl
       << fseq << endl;
}

void MRFShowProcessor::show_prof(ostream& os, const MRF& model) {
    Eigen::IOFormat fmt(prec, Eigen::DontAlignCols, sep, sep);
    size_t n = model.get_length();
    string seq = model.get_seq();
    os << "# PROFILE" << endl
       << "RES"     << sep << format_symbol() << endl;
    for (size_t i = 0; i < n; ++i)
        os << seq[i] << " " << i + 1 << sep << model.get_psfm(i).format(fmt) << endl;
}

void MRFShowProcessor::show_mrf(ostream& os, const MRF& model) {
    Eigen::IOFormat fmt(prec, Eigen::DontAlignCols, sep, sep);
    size_t n = model.get_length();
    string seq = model.get_seq();
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    os << "# NODE"  << endl
       << "RES"     << sep << format_symbol() << endl;
    for (size_t i = 0; i < n; ++i)
        os << seq[i] << " " << i + 1 << sep << model.get_node(i).get_weight().format(fmt) << endl;
    os << "# EDGE"  << endl
       << "RES1"    << sep << "RES2" << sep << format_paired_symbol() << endl;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        size_t i = pos->idx1;
        size_t j = pos->idx2;
        os << seq[i] << " " << i + 1 << sep
           << seq[j] << " " << j + 1 << sep
           << model.get_edge(i, j).get_weight().format(fmt) << endl;
    }
}

void MRFShowProcessor::show_mrf_node(ostream& os, const MRF& model, const size_t& pos) {
    const string symbol = abc.get_canonical();
    VectorXf v = model.get_node(pos).get_weight();
    os << setprecision(prec);
    for (size_t i = 0; i < symbol.size(); ++i)
        os << symbol[i] << sep << v(i) << endl;
}

void MRFShowProcessor::show_mrf_edge(ostream& os, const MRF& model, const EdgeIndex& pos) {
    Eigen::IOFormat fmt(prec, 0, "\t");
    os << setprecision(prec);
    const string symbol = abc.get_canonical();
    if (model.has_edge(pos)) {
        MatrixXf w = model.get_edge(pos).get_weight();
        if (pos.idx1 > pos.idx2) w.transposeInPlace();
        os << pos.idx1 << "\\" << pos.idx2;
        for (size_t i = 0; i < symbol.size(); ++i) os << sep << symbol[i];
        os << endl;
        for (size_t i = 0; i < symbol.size(); ++i)
            os << symbol[i] << sep << w.row(i).format(fmt) << endl;
    } else {
        os << "Edge (" << pos.idx1 << ", " << pos.idx2 << ") does not exist" << endl;
    }
}

string MRFShowProcessor::format_symbol() {
    const string symbol = abc.get_canonical();
    string ret;
    ret += symbol[0];
    for (size_t i = 1; i < symbol.size(); ++i) ret += sep + symbol[i];
    return ret;
}

string MRFShowProcessor::format_paired_symbol() {
    const string symbol = abc.get_canonical();
    string ret;
    ret += symbol[0];
    ret += symbol[0];
    for (size_t i = 0; i < symbol.size(); ++i)
        for (size_t j = 0; j < symbol.size(); ++j)
            if (i != 0 || j != 0) ret += sep + symbol[i] + symbol[j];
    return ret;
}

ostream& operator<<(ostream& os, const MRF& model) {
    const string sep = "\t";
    os << "# PMRF INFO" << endl
       << "FMTNUM"  << sep << model.get_fmtnum() << endl
       << "NEFF"    << sep << model.get_neff() << endl
       << "LENGTH"  << sep << model.get_length() << endl
       << "EDGES"   << sep << model.get_num_edge() << endl;
   return os;
}
