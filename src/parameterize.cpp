#include "parameterize.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <memory>
#include <limits>

#include "util/common.h"
#include "seq/profile.h"
#include "seq/bgfreq.h"
#include "seq/subsmat.h"

#include "core.h"

using std::make_shared;

/**
   @class MRFParameterizer::Parameter
 */

void MRFParameterizer::Parameter::init(const MRF& model, const Option& opt) {
    num_var = model.get_num_var();
    n_node = length * num_var;
    num_edge = eidx.size();
    num_var2 = num_var * num_var;
    n_edge = num_edge * num_var2;
    n = n_node + n_edge;
    LBFGS::Parameter::malloc_x();       // This should be called after setting this->n
    set_opt(opt);
}

void MRFParameterizer::Parameter::set_opt(const Option& opt) {
    /* L-BFGS options */
    LBFGS::Parameter::init_param();     // This should be called before setting parameters
    opt_param.m = opt.corr;                                 // number of correction
    opt_param.epsilon = (lbfgsfloatval_t) opt.epsilon;      // convergence criterion
    opt_param.past = opt.past;                              // stopping criterion
    opt_param.delta = (lbfgsfloatval_t) opt.delta;          // stopping criterion
    opt_param.max_iterations = opt.max_iterations;          // maximum iteration
    opt_param.linesearch = opt.linesearch;
    //opt_param.max_linesearch = 50;
}


/**
   @class MRFParameterizer::SymmParameter
 */

MRFParameterizer::SymmParameter::SymmParameter(const MRF& model, const Option& opt) : Parameter(model, opt) {
    edge_idxs = model.get_edge_idxs();
    eidx.clear();
    int i = -1;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos)
        eidx[EdgeIndex(pos->idx1, pos->idx2)] = ++i;
    init(model, opt);
}


/**
   @class MRFParameterizer::LocalParameter
 */

MRFParameterizer::LocalParameter::LocalParameter(const MRF& model, const Option& opt, const Parameter& p, const int& obs_node) : Parameter(model, opt), obs_node(obs_node) {
    eidx.clear();
    int ei = -1;
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = p.eidx.cbegin(); pos != p.eidx.cend(); ++pos)
        if (pos->first.idx1 == obs_node) eidx[EdgeIndex(pos->first.idx1, pos->first.idx2)] = ++ei;
    num_var = model.get_num_var();
    n_node = 1 * num_var;
    num_edge = eidx.size();
    num_var2 = num_var * num_var;
    n_edge = num_edge * num_var2;
    n = n_node + n_edge;
    LBFGS::Parameter::malloc_x();       // This should be called after setting this->n
    set_opt(opt);
    // set x values
    int i = p.get_nidx(obs_node, 0);
    for (int k = 0; k < num_var; ++k) x[k] = p.x[i + k];
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = eidx.cbegin(); pos != eidx.cend(); ++pos) {
        i = get_eidx_edge(pos->second, 0, 0);
        int j = p.get_eidx(pos->first.idx1, pos->first.idx2, 0, 0);
        for (int k = 0; k < num_var2; ++k) x[i + k] = p.x[j + k];
    }
}

/**
   @class MRFParameterizer::AsymParameter
 */

MRFParameterizer::AsymParameter::AsymParameter(const MRF& model, const Option& opt) : Parameter(model, opt) {
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    eidx.clear();
    int i = -1;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        eidx[EdgeIndex(pos->idx1, pos->idx2)] = ++i;
        eidx[EdgeIndex(pos->idx2, pos->idx1)] = ++i;
    }
    init(model, opt);
}

void MRFParameterizer::AsymParameter::set_x(const LocalParameter& p, const int& obs_node) {
    int i, j;
    i = get_nidx(obs_node, 0);
    j = p.get_nidx(obs_node, 0);
    for (int k = 0; k < num_var; ++k) x[i + k] = p.x[j + k];
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = p.eidx.cbegin(); pos != p.eidx.cend(); ++pos) {
        i = get_eidx(pos->first.idx1, pos->first.idx2, 0, 0);
        j = p.get_eidx(pos->first.idx1, pos->first.idx2, 0, 0);
        for (int k = 0; k < num_var2; ++k) x[i + k] = p.x[j + k];
    }
}


/**
   @class MRFParameterizer::L2Regularization
 */

lbfgsfloatval_t MRFParameterizer::L2Regularization::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    regularize_node(x, g, fx);
    regularize_edge(x, g, fx);
    return fx;
}

void MRFParameterizer::L2Regularization::regularize_node(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const float& lambda = opt.lambda1;
    for (int k = param.nidx_beg(); k < param.nidx_end(); ++k) {
        fx += lambda * square(x[k]);
        g[k] += 2. * lambda * x[k];
    }
}

void MRFParameterizer::L2Regularization::regularize_edge(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const float& lambda = opt.lambda2;
    const float twolambda = 2. * lambda;
    for (int k = param.eidx_beg(); k < param.eidx_end(); ++k) {
        fx += lambda * square(x[k]);
        g[k] += twolambda * x[k];
    }
}


/**
   @class MRFParameterizer::Pseudolikelihood
 */

MRFParameterizer::Pseudolikelihood::Pseudolikelihood(const TraceVector& traces, const SymmParameter& param, const VectorXf& seq_weight, const float& neff)
: traces(traces), 
  param(param), 
  seq_weight(seq_weight),
  neff(neff) {
    data.resize(traces.size(), param.length);
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        for (size_t j = 0; j < param.length; ++j)
            data(i, j) = param.abc.get_idx(seq[j]);
    }
}

lbfgsfloatval_t MRFParameterizer::Pseudolikelihood::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        const float sw = seq_weight(i);
        const MatrixXf logpot = calc_logpot(x, i, seq);
        const VectorXf logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, i, seq, sw);
        update_gradient(x, g, logpot, logz, i, seq, sw);
    }
    return fx;
}

MatrixXf MRFParameterizer::Pseudolikelihood::calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq) {
    const Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > nw(x + param.nidx_beg(), param.length, param.num_var);
    MatrixXf logpot(nw.transpose());
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.cbegin(); pos != param.eidx.cend(); ++pos) {
        const int i = pos->first.idx1;
        const int j = pos->first.idx2;
        int ei = pos->second;
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

VectorXf MRFParameterizer::Pseudolikelihood::logsumexp(const MatrixXf& b) {
    VectorXf B = b.rowwise().maxCoeff();
    VectorXf r = VectorXf::Zero(b.rows());
    for (int j = 0; j < b.cols(); ++j)
        r += (b.col(j) - B).unaryExpr(&exp);
    return r.unaryExpr(&log) + B;
}

VectorXf MRFParameterizer::Pseudolikelihood::calc_logz(const MatrixXf& logpot) {
    return logsumexp(logpot.transpose());
}

void MRFParameterizer::Pseudolikelihood::update_obj_score(lbfgsfloatval_t& fx, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw) {
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

void MRFParameterizer::Pseudolikelihood::update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXf& logpot, const VectorXf& logz, const size_t m, const string& seq, const float& sw) {
    string letters = param.abc.get_canonical();
    float w;
    MatrixXf nodebel(param.num_var, param.length);
    for (int t = 0; t < param.num_var; ++t)
        nodebel.row(t) = (logpot.row(t) - logz.transpose()).unaryExpr(&exp);
    nodebel *= sw;
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > dv(g, param.length, param.num_var);
    for (int i = 0; i < param.length; ++i) {
        int xi = data(m, i);
        if (xi < param.num_var) {
            dv(i, xi) -= sw;
        } else {
            string s = param.abc.get_degeneracy(seq[i], &w);
            for (string::iterator c = s.begin(); c != s.end(); ++c) {
                xi = param.abc.get_idx(*c);
                dv(i, xi) -= w * sw;
            }
        }
        dv.row(i) += nodebel.col(i);
    }
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.cbegin(); pos != param.eidx.cend(); ++pos) {
        const int i = pos->first.idx1;
        const int j = pos->first.idx2;
        const int ei = pos->second;
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
   @class MRFParameterizer::AsymPseudolikelihood
 */

MRFParameterizer::AsymPseudolikelihood::AsymPseudolikelihood(const TraceVector& traces, const Parameter& param, const VectorXf& seq_weight, const float& neff)
: traces(traces), 
  param(param), 
  lp(NULL),
  seq_weight(seq_weight),
  neff(neff) {
    data.resize(traces.size(), param.length);
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        for (size_t j = 0; j < param.length; ++j)
            data(i, j) = param.abc.get_idx(seq[j]);
    }
}

lbfgsfloatval_t MRFParameterizer::AsymPseudolikelihood::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    const Parameter* p;
    if (lp != NULL) p = lp;
    else p = &param;
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        const float sw = seq_weight(i);
        const VectorXf logpot = calc_logpot(x, i, seq, *p);
        const float logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, i, seq, sw);
        update_gradient(x, g, logpot, logz, i, seq, sw, *p);
    }
    return fx;
}

VectorXf MRFParameterizer::AsymPseudolikelihood::calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq, const Parameter& param) {
    const Eigen::Map<const VectorXf> nw(x + param.get_nidx(obs_node, 0), param.num_var);
    VectorXf logpot(nw);
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.cbegin(); pos != param.eidx.cend(); ++pos) {
        if (pos->first.idx1 != obs_node) continue;
        const int i = pos->first.idx1;
        const int j = pos->first.idx2;
        const int xi = data(m, i);
        const int xj = data(m, j);
        const Eigen::Map<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > ew(x + param.get_eidx_edge(pos->second, 0, 0), param.num_var, param.num_var);
        if (xj < param.num_var) {
            logpot += ew.col(xj);
        } else {
            float wj;
            string sj = param.abc.get_degeneracy(seq[j], &wj);
            for (string::iterator c = sj.begin(); c != sj.end(); ++c)
                logpot += wj * ew.col(param.abc.get_idx(*c));
        }
    }
    return logpot;
}

float MRFParameterizer::AsymPseudolikelihood::calc_logz(const VectorXf& logpot) {
    return log(logpot.unaryExpr(&exp).sum());
}

void MRFParameterizer::AsymPseudolikelihood::update_obj_score(lbfgsfloatval_t& fx, const VectorXf& logpot, const float& logz, const size_t m, const string& seq, const float& sw) {
    int i = obs_node;
    const int xi = data(m, i);
    if (xi < param.num_var) {
        fx += -sw * (logpot(xi) - logz);
    } else {
        FloatType w;
        string s = param.abc.get_degeneracy(seq[i], &w);
        for (string::iterator c = s.begin(); c != s.end(); ++c)
            fx += w * (-sw * logpot(param.abc.get_idx(*c)) + sw * logz);
    }
}

void MRFParameterizer::AsymPseudolikelihood::update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const VectorXf& logpot, const float& logz, const size_t m, const string& seq, const float& sw, const Parameter& param) {
    string letters = param.abc.get_canonical();
    float w;
    VectorXf nodebel(param.num_var);
    nodebel = (logpot - VectorXf::Constant(param.num_var, logz)).unaryExpr(&exp);
    nodebel *= sw;
    Eigen::Map<VectorXf> dv(g + param.get_nidx(obs_node, 0), param.num_var);
    int i = obs_node;
    int xi = data(m, i);
    if (xi < param.num_var) {
        dv(xi) -= sw;
    } else {
        string s = param.abc.get_degeneracy(seq[i], &w);
        for (string::iterator c = s.begin(); c != s.end(); ++c) {
            xi = param.abc.get_idx(*c);
            dv(xi) -= w * sw;
        }
    }
    dv += nodebel;
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.cbegin(); pos != param.eidx.cend(); ++pos) {
        if (pos->first.idx1 != obs_node) continue;
        const int i = pos->first.idx1;
        const int j = pos->first.idx2;
        int xi = data(m, i);
        int xj = data(m, j);
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > dw(g + param.get_eidx_edge(pos->second, 0, 0), param.num_var, param.num_var);
        if (xi < param.num_var && xj < param.num_var) {
            dw(xi, xj) -= sw;
            dw.col(xj) += nodebel;
        } else {
            float wi, wj;
            string si = param.abc.get_degeneracy(seq[i], &wi);
            string sj = param.abc.get_degeneracy(seq[j], &wj);
            w = wi * wj;
            for (string::iterator ci = si.begin(); ci != si.end(); ++ci) {
                xi = param.abc.get_idx(*ci);
                for (string::iterator cj = sj.begin(); cj != sj.end(); ++cj) {
                    xj = param.abc.get_idx(*cj);
                    dw(xi, xj) -= w * sw;
                    for (int t = 0; t < param.num_var; ++t) {
                        dw.col(xj) += w * nodebel;
                    }
                }
            }
        }
    }
}

/**
   @class MRFParameterizer::ObjectiveFunction
 */

lbfgsfloatval_t MRFParameterizer::ObjectiveFunction::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (int i = 0; i < param.n; ++i) g[i] = 0.;
    for (vector<PtrObjFunc>::const_iterator pos = funcs.cbegin(); pos != funcs.cend(); ++pos)
        fx += (*pos)->evaluate(x, g, n, step);
    return fx;
}

/**
   @class MRFParameterizer
 */

int MRFParameterizer::parameterize(MRF& model, const TraceVector& traces) {
    size_t length = model.get_length();
    size_t num_var = model.get_num_var();
    size_t num_edge = model.get_num_edge();
    int ret = 0;
    /* sequeing weights */
    vector<string> msa = traces.get_matched_aseq_vec();
    msa = msa_analyzer.termi_gap_remover->filter(msa);
    VectorXf sw = msa_analyzer.seq_weight_estimator->estimate(msa);
    sw /= sw.sum();
    /* effective number of sequences */
    float neff = msa_analyzer.eff_seq_num_estimator->estimate(msa);
    sw *= neff;
    /* sequence profile */
    MatrixXf psfm = MatrixXf::Zero(length, num_var);
    FloatType gap_prob = 0.;
    psfm.leftCols<20>() = calc_profile(traces) * (1. - gap_prob);
    psfm.rightCols<1>().setConstant(gap_prob);
    model.set_psfm(psfm);
    /* MRF */
    if (opt.regedge_sc_deg) opt.regedge_lambda *= 2. * (float) num_edge / (float) length;
    if (opt.regedge_sc_neff) opt.regedge_lambda *= neff;
    if (opt.asymmetric) opt.regedge_lambda *= 0.5;
    if (opt.regul == RegulMethod::REGUL_L2) {
        opt.l2_opt.lambda1 = opt.regnode_lambda;
        opt.l2_opt.lambda2 = opt.regedge_lambda;
    }
    vector<PtrObjFunc> funcs;
    LBFGS::Optimizer optimizer;
    if (opt.asymmetric) {
        AsymParameter param(model, optim_opt);
        AsymPseudolikelihood pll(traces, param, sw, neff);
        for (size_t i = 0; i < length; ++i) {
            //std::cout << "obs_node = " << i << std::endl;
            LocalParameter local_param(model, optim_opt, param, i);
            pll.set_local_param(&local_param);
            pll.set_obs_node(i);
            funcs.push_back(&pll);
            L2Regularization l2(local_param, opt.l2_opt);
            if (opt.regul == RegulMethod::REGUL_L2) funcs.push_back(&l2);
            ObjectiveFunction obj_func(funcs, local_param);
            int r = optimizer.optimize(&local_param, &obj_func);
            funcs.clear();
            if (r < ret) ret = r;
            param.set_x(local_param, i);
        }
        adjust_gauge(param);
        update_model(model, param);
    } else {
        SymmParameter param(model, optim_opt);
        Pseudolikelihood pll(traces, param, sw, neff);
        funcs.push_back(&pll);
        L2Regularization l2(param, opt.l2_opt);
        if (opt.regul == RegulMethod::REGUL_L2) funcs.push_back(&l2);
        ObjectiveFunction obj_func(funcs, param);
        ret = optimizer.optimize(&param, &obj_func);
        adjust_gauge(param);
        update_model(model, param);
    }
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

void MRFParameterizer::update_model(MRF& model, const SymmParameter& param) {
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

void MRFParameterizer::update_model(MRF& model, const AsymParameter& param) {
    string letters = param.abc.get_canonical();
    for (size_t i = 0; i < model.get_length(); ++i) {
        VectorXf w = VectorXf::Zero(param.num_var);
        for (int t = 0; t < param.num_var; ++t) w(t) = (FloatType)param.x[param.get_nidx(i, letters[t])];
        model.get_node(i).set_weight(w);
    }
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    for (EdgeIndexVector::const_iterator pos = edge_idxs.cbegin(); pos != edge_idxs.cend(); ++pos) {
        const int i = pos->idx1;
        const int j = pos->idx2;
        MatrixXf w = MatrixXf::Zero(param.num_var, param.num_var);
        for (int p = 0; p < param.num_var; ++p)
            for (int q = 0; q < param.num_var; ++q)
                w(p, q) += 0.5 * (param.x[param.get_eidx(i, j, p, q)] + param.x[param.get_eidx(j, i, q, p)]);
        model.get_edge(i, j).set_weight(w);
    }
}

void MRFParameterizer::adjust_gauge(SymmParameter& param) {
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > v(param.x, param.length, param.num_var);
    for (size_t i = 0; i < param.length; ++i) {
        v.row(i) -= VectorXf::Constant(param.num_var, v.mean());
    }
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.cbegin(); pos != param.eidx.cend(); ++pos) {
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > w(param.x + param.get_eidx_edge(pos->second, 0, 0), param.num_var, param.num_var);
        float totmn = w.mean();
        w += MatrixXf::Constant(param.num_var, param.num_var, totmn);
    }
}

void MRFParameterizer::adjust_gauge(AsymParameter& param) {
    Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > v(param.x, param.length, param.num_var);
    for (size_t i = 0; i < param.length; ++i) {
        v.row(i) -= VectorXf::Constant(param.num_var, v.mean());
    }
    for (std::unordered_map<EdgeIndex, int>::const_iterator pos = param.eidx.cbegin(); pos != param.eidx.cend(); ++pos) {
        const int i = pos->first.idx1;
        Eigen::Map<Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> > w(param.x + param.get_eidx_edge(pos->second, 0, 0), param.num_var, param.num_var);
        VectorXf rowmn(w.rowwise().mean());
        VectorXf colmn(w.colwise().mean());
        float totmn = w.mean();
        v.row(i) += rowmn;
        v.row(i) -= VectorXf::Constant(param.num_var, totmn);
        for (size_t k = 0; k < param.num_var; ++k) {
            w.col(k) -= rowmn;
            w.row(k) -= colmn;
        }
        w += MatrixXf::Constant(param.num_var, param.num_var, totmn);
    }
}
