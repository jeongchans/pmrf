#include "parameterize.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "util/common.h"
#include "seq/profile.h"
#include "seq/bgfreq.h"
#include "seq/subsmat.h"

#include "core.h"

/**
   @class MRFParameterizer::Parameter
 */

MRFParameterizer::Parameter::Parameter(const MRF& model, const Option& opt) : abc(model.get_alphabet()), length(model.get_length()) {
    edge_idxs = model.get_edge_idxs();
    eidx.clear();
    int i = 0;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        eidx[EdgeIndex(pos->idx1, pos->idx2)] = i;
        i++;
    }
    num_var = model.get_num_var();
    n_node = length * num_var;
    num_edge = eidx.size();
    num_var2 = num_var * num_var;
    n_edge = num_edge * num_var2;
    n = n_node + n_edge;
    init();
    opt_param.linesearch = opt.linesearch;
    opt_param.delta = (lbfgsfloatval_t)(opt.delta);
    opt_param.past = opt.past;
    opt_param.max_iterations = opt.max_iterations;
}

int MRFParameterizer::Parameter::get_nidx(const int& i, const char& p) const {
    return i * num_var + abc.get_idx(p);
}

inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const char& q) const {
    return get_eidx(i, j, abc.get_idx(p), abc.get_idx(q));
}

inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const int& xi, const char& q) const {
    return get_eidx(i, j, xi, abc.get_idx(q));
}

inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const int& xj) const {
    return get_eidx(i, j, abc.get_idx(p), xj);
}

inline int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const int& xi, const int& xj) const {
    return get_eidx_edge(eidx.at(EdgeIndex(i, j)), xi, xj);
}

inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const char& p, const char& q) const {
    return get_eidx_edge(ei, abc.get_idx(p), abc.get_idx(q));
}

inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const int& xi, const char& q) const {
    return get_eidx_edge(ei, xi, abc.get_idx(q));
}

inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const char& p, const int& xj) const {
    return get_eidx_edge(ei, abc.get_idx(p), xj);
}

inline int MRFParameterizer::Parameter::get_eidx_edge(const int& ei, const int& xi, const int& xj) const {
    return eidx_beg() + ei * num_var2 + xi * num_var + xj;
}

/**
   @class MRFParameterizer::L2Regularization
 */

void MRFParameterizer::L2Regularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    regularize_node(x, g, fx);
    regularize_edge(x, g, fx);
}

void MRFParameterizer::L2Regularization::regularize_node(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const float& lambda = opt.lambda1;
    for (int k = param.nidx_beg(); k < param.nidx_end(); ++k) {
        fx += lambda * square(x[k]);
        g[k] += 2. * lambda * x[k];
    }
}

void MRFParameterizer::L2Regularization::regularize_edge(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    float lambda = opt.lambda2;
    if (opt.sc) lambda *= 2. * ((float)(param.eidx.size())) / ((float)(param.length));
    float twolambda = 2. * lambda;
    for (int k = param.eidx_beg(); k < param.eidx_end(); ++k) {
        fx += lambda * square(x[k]);
        g[k] += twolambda * x[k];
    }
}

/**
   @class MRFParameterizer::ObjectiveFunction
 */

MRFParameterizer::ObjectiveFunction::ObjectiveFunction(const TraceVector& traces, Parameter& param, Option& opt, const MSAAnalyzer& msa_analyzer)
: traces(traces), 
  param(param), 
  opt(opt), 
  msa_analyzer(msa_analyzer),
  l2_func(param, opt.l2_opt) {
    data.resize(traces.size(), param.length);
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        for (size_t j = 0; j < param.length; ++j)
            data(i, j) = param.abc.get_idx(seq[j]);
    }
    vector<string> msa = traces.get_matched_aseq_vec();
    msa = msa_analyzer.termi_gap_remover->filter(msa);
    seq_weight.resize(traces.size());
    seq_weight = msa_analyzer.seq_weight_estimator->estimate(msa);
    seq_weight /= seq_weight.sum();
    seq_weight *= msa_analyzer.eff_seq_num_estimator->estimate(msa);
}

lbfgsfloatval_t MRFParameterizer::ObjectiveFunction::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (int i = 0; i < param.n; ++i) g[i] = 0.;
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        const float sw = seq_weight(i);
        const MatrixXf logpot = calc_logpot(x, i, seq);
        const VectorXf logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, i, seq, sw);
        update_gradient(x, g, logpot, logz, i, seq, sw);
    }
    if (opt.regul == RegulMethod::RegulMethod::L2) l2_func.regularize(x, g, fx);
    return fx;
}

MatrixXf MRFParameterizer::ObjectiveFunction::calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq) {
    const Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > nw(x + param.nidx_beg(), param.length, param.num_var);
    MatrixXf logpot(nw.transpose());
    for (int ei = 0; ei < param.num_edge; ++ei) {
        const int i = param.edge_idxs[ei].idx1;
        const int j = param.edge_idxs[ei].idx2;
        const int xi = data(m, i);
        const int xj = data(m, j);
        const Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ew(x + param.get_eidx_edge(ei, 0, 0), param.num_var, param.num_var);
        if (xi < param.num_var) {
            logpot.col(j) += ew.row(xi);
        } else {
            float wi;
            string si = param.abc.get_degeneracy(seq[i], &wi);
            for (string::iterator c = si.begin(); c != si.end(); ++c)
                logpot.col(j) += wi * ew.row(param.abc.get_idx(*c));
        }
        if (xj < param.num_var) {
            logpot.col(i) += ew.col(xj);
        } else {
            float wj;
            string sj = param.abc.get_degeneracy(seq[j], &wj);
            for (string::iterator c = sj.begin(); c != sj.end(); ++c)
                logpot.col(i) += wj * ew.col(param.abc.get_idx(*c));
        }
    }
    return logpot;
}

VectorXf MRFParameterizer::ObjectiveFunction::logsumexp(const MatrixXf& b) {
    VectorXf B = b.rowwise().maxCoeff();
    VectorXf r = VectorXf::Zero(b.rows());
    for (int j = 0; j < b.cols(); ++j)
        r += (b.col(j) - B).unaryExpr(&exp);
    return r.unaryExpr(&log) + B;
}

VectorXf MRFParameterizer::ObjectiveFunction::calc_logz(const MatrixXf& logpot) {
    return logsumexp(logpot.transpose());
}

void MRFParameterizer::ObjectiveFunction::update_obj_score(lbfgsfloatval_t& fx, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw) {
    for (int i = 0; i < param.length; ++i) {
        const int xi = data(m, i);
        if (xi < param.num_var) {
            fx += -sw * (logpot(xi, i) - logz(i));
        } else {
            FloatType w;
            string s = param.abc.get_degeneracy(seq[i], &w);
            for (string::iterator c = s.begin(); c != s.end(); ++c)
                fx += w * (-sw * logpot(param.abc.get_idx(*c), i) + sw * logz(i));
        }
    }
}

void MRFParameterizer::ObjectiveFunction::update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw) {
    string letters = param.abc.get_canonical();
    float w;
    MatrixXf nodebel(param.num_var, param.length);
    for (int t = 0; t < param.num_var; ++t)
        nodebel.row(t) = (logpot.row(t) - logz.transpose()).unaryExpr(&exp);
    nodebel *= sw;
    for (int i = 0; i < param.length; ++i) {
        string s = param.abc.get_degeneracy(seq[i], &w);
        for (string::iterator c = s.begin(); c != s.end(); ++c)
            g[param.get_nidx(i, *c)] -= w * sw;
        for (int t = 0; t < param.num_var; ++t)
            g[param.get_nidx(i, letters[t])] += nodebel(t, i);
    }
    for (int ei = 0; ei < param.num_edge; ++ei) {
        const int i = param.edge_idxs[ei].idx1;
        const int j = param.edge_idxs[ei].idx2;
        int xi = data(m, i);
        int xj = data(m, j);
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > dw(g + param.get_eidx_edge(ei, 0, 0), param.num_var, param.num_var);
        if (xi < param.num_var && xj < param.num_var) {
            dw(xi, xj) -= sw * 2.;
            dw.col(xj) += nodebel.col(i);
            dw.row(xi) += nodebel.col(j);
        } else {
            float wi, wj;
            string si = param.abc.get_degeneracy(seq[i], &wi);
            string sj = param.abc.get_degeneracy(seq[j], &wj);
            w = wi * wj;
            for (string::iterator ci = si.begin(); ci != si.end(); ++ci) {
                xi = param.abc.get_idx(*ci);
                for (string::iterator cj = sj.begin(); cj != sj.end(); ++cj) {
                    xj = param.abc.get_idx(*cj);
                    dw(xi, xj) -= w * sw * 2.;
                    for (int t = 0; t < param.num_var; ++t) {
                        dw.col(xj) += w * nodebel.col(i);
                        dw.row(xi) += w * nodebel.col(j);
                    }
                }
            }
        }
    }
}

/**
   @class MRFParameterizer
 */

int MRFParameterizer::parameterize(MRF& model, const TraceVector& traces) {
    Parameter param(model, optim_opt);
    MatrixXf psfm = MatrixXf::Zero(param.length, param.num_var);
    FloatType gap_prob = 0.;
    psfm.leftCols<20>() = calc_profile(traces) * (1. - gap_prob);
    psfm.rightCols<1>().setConstant(gap_prob);
    ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);
    LBFGS::Optimizer optimizer;
    int ret = optimizer.optimize(&param, &obj_func);
    model.set_psfm(psfm);
    update_model(model, param);
    return ret;
}

MatrixXf MRFParameterizer::calc_profile(const TraceVector& traces) {
    AminoAcid aa;
    ProfileBuilder profile_builder(aa);
    SMMEmitProbEstimator emit_prob_estimator(ROBINSON_BGFREQ.get_array(aa), 
                                             BLOSUM62_MATRIX.get_array(aa), 
                                             3.2);
    profile_builder.set_emit_prob_estimator(&emit_prob_estimator);
    profile_builder.set_seq_weight_estimator(msa_analyzer.seq_weight_estimator.get());
    profile_builder.set_eff_seq_num_estimator(msa_analyzer.eff_seq_num_estimator.get());
    profile_builder.set_msa_filter(msa_analyzer.termi_gap_remover.get());
    Profile profile = profile_builder.build(traces);
    return profile.get_prob();
}

void MRFParameterizer::update_model(MRF& model, Parameter& param) {
    string letters = param.abc.get_canonical();
    for (size_t i = 0; i < model.get_length(); ++i) {
        VectorXf w = VectorXf::Zero(param.num_var);
        for (int t = 0; t < param.num_var; ++t) w(t) = (FloatType)param.x[param.get_nidx(i, letters[t])];
        model.get_node(i).set_weight(w);
    }
    for (unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.begin(); pos != param.eidx.end(); ++pos) {
        int i = pos->first.idx1;
        int j = pos->first.idx2;
        MatrixXf w = MatrixXf::Zero(param.num_var, param.num_var);
        for (int p = 0; p < param.num_var; ++p)
            for (int q = 0; q < param.num_var; ++q) w(p, q) = (FloatType)param.x[param.get_eidx(i, j, letters[p], letters[q])];
        model.get_edge(i, j).set_weight(w);
    }
}
