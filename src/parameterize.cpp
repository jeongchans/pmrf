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

MatrixXi trace_to_var(const TraceVector& traces, const Alphabet& abc) {
    const size_t row = traces.size();
    const size_t col = traces[0].get_matched_aseq().size();
    MatrixXi data(row, col);
    for (size_t i = 0; i < row; ++i) {
        string seq = traces[i].get_matched_aseq();
        for (size_t j = 0; j < col; ++j)
            data(i, j) = abc.get_idx(seq[j]);
    }
    return data;
}

/* log-sum-exp trick (https://en.wikipedia.org/wiki/LogSumExp) */
lbfgsfloatval_t logsumexp(const VectorXl& x) {
    lbfgsfloatval_t b = x.maxCoeff();
    return log((x - VectorXl::Constant(x.size(), b)).unaryExpr(&exp).sum()) + b;
}

/**
   @class MRFParameterizer::Parameter
 */

void MRFParameterizer::Parameter::init(const vector<NodeIndex>& v_idxs, const vector<EdgeIndex>& w_idxs, const Option& opt) {
    num_node = v_idxs.size();
    n_v = num_node * num_var;
    num_edge = w_idxs.size();
    n_w = num_edge * num_var2;
    n = n_v + n_w;
    size_t offset = 0;
    v_offset.clear();
    for (auto it = v_idxs.cbegin(); it != v_idxs.cend(); ++it, offset += num_var) v_offset[*it] = offset;
    w_offset.clear();
    for (auto it = w_idxs.begin(); it != w_idxs.end(); ++it, offset += num_var2) w_offset[*it] = offset;
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
    opt_param.max_linesearch = opt.max_linesearch;
}


/**
   @class MRFParameterizer::SymmParameter
 */

MRFParameterizer::SymmParameter::SymmParameter(const MRF& model, const Option& opt) 
: Parameter(model.get_alphabet(), model.get_num_var()) {
    vector<NodeIndex> v_idxs;
    for (size_t i = 0; i < model.get_length(); ++i) v_idxs.push_back(i);
    vector<EdgeIndex> w_idxs = model.get_edge_idxs();
    init(v_idxs, w_idxs, opt);
}


/**
   @class MRFParameterizer::LocalParameter
 */

MRFParameterizer::LocalParameter::LocalParameter(const MRF& model, const Option& opt, const Parameter& p, const size_t& obs_node) 
: Parameter(model.get_alphabet(), model.get_num_var()) {
    vector<NodeIndex> v_idxs(1, obs_node);
    vector<EdgeIndex> w_idxs;
    for (auto it = p.w_offset.cbegin(); it != p.w_offset.cend(); ++it)
        if (it->first.idx1 == obs_node) w_idxs.push_back(it->first);
    init(v_idxs, w_idxs, opt);
    // set x values
    size_t offset1 = v_offset.at(obs_node);
    size_t offset2 = p.v_offset.at(obs_node);
    for (size_t k = 0; k < num_var; ++k, ++offset1, ++offset2) x[offset1] = p.x[offset2];
    for (auto it = w_offset.cbegin(); it != w_offset.cend(); ++it) {
        offset1 = it->second;
        offset2 = p.w_offset.at(it->first);
        for (size_t k = 0; k < num_var2; ++k, ++offset1, ++offset2) x[offset1] = p.x[offset2];
    }
}


/**
   @class MRFParameterizer::AsymParameter
 */

MRFParameterizer::AsymParameter::AsymParameter(const MRF& model, const Option& opt)
: Parameter(model.get_alphabet(), model.get_num_var()) {
    vector<NodeIndex> v_idxs;
    for (size_t i = 0; i < model.get_length(); ++i) v_idxs.push_back(i);
    vector<EdgeIndex> w_idxs;
    vector<EdgeIndex> sym_w_idxs = model.get_edge_idxs();
    for (auto it = sym_w_idxs.cbegin(); it != sym_w_idxs.cend(); ++it) {
        w_idxs.push_back(EdgeIndex(it->idx1, it->idx2));
        w_idxs.push_back(EdgeIndex(it->idx2, it->idx1));
    }
    init(v_idxs, w_idxs, opt);
}

void MRFParameterizer::AsymParameter::set_x(const LocalParameter& p) {
    for (auto it = p.v_offset.cbegin(); it != p.v_offset.cend(); ++it) {
        size_t offset1 = v_offset[it->first];
        size_t offset2 = it->second;
        for (size_t k = 0; k < num_var; ++k) { x[offset1] = p.x[offset2]; ++offset1; ++offset2; }
    }
    for (auto it = p.w_offset.cbegin(); it != p.w_offset.cend(); ++it) {
        size_t offset1 = w_offset[it->first];
        size_t offset2 = it->second;
        for (size_t k = 0; k < num_var2; ++k) { x[offset1] = p.x[offset2]; ++offset1; ++offset2; }
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
    const lbfgsfloatval_t& lambda = opt.lambda1;
    const lbfgsfloatval_t twolambda = 2. * lambda;
    for (int k = param.v_beg_pos(); k < param.v_end_pos(); ++k) {
        fx += lambda * square(x[k]);
        g[k] += twolambda * x[k];
    }
}

void MRFParameterizer::L2Regularization::regularize_edge(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const lbfgsfloatval_t& lambda = opt.lambda2;
    const lbfgsfloatval_t twolambda = 2. * lambda;
    for (int k = param.w_beg_pos(); k < param.w_end_pos(); ++k) {
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
    data = trace_to_var(traces, param.abc);
}

lbfgsfloatval_t MRFParameterizer::Pseudolikelihood::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        const float sw = seq_weight(i);
        const MatrixXl logpot = calc_logpot(x, i, seq);
        const VectorXl logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, i, seq, sw);
        update_gradient(x, g, logpot, logz, i, seq, sw);
    }
    return fx;
}

MatrixXl MRFParameterizer::Pseudolikelihood::calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq) {
    MapKMatrixXl v(x + param.v_beg_pos(), param.num_node, param.num_var);
    MatrixXl logpot(v.transpose());
    for (auto it = param.w_offset.cbegin(); it != param.w_offset.cend(); ++it) {
        const int i = it->first.idx1;
        const int j = it->first.idx2;
        MapKMatrixXl w(x + it->second, param.num_var, param.num_var);
        const int xi = data(m, i);
        const int xj = data(m, j);
        if (xi < param.num_var) logpot.col(j) += w.row(xi);
        else {
            float wi;
            string si = param.abc.get_degeneracy(seq[i], &wi);
            for (auto c = si.cbegin(); c != si.cend(); ++c) logpot.col(j) += wi * w.row(param.abc.get_idx(*c));
        }
        if (xj < param.num_var) logpot.col(i) += w.col(xj);
        else {
            float wj;
            string sj = param.abc.get_degeneracy(seq[j], &wj);
            for (auto c = sj.cbegin(); c != sj.cend(); ++c) logpot.col(i) += wj * w.col(param.abc.get_idx(*c));
        }
    }
    return logpot;
}

VectorXl MRFParameterizer::Pseudolikelihood::calc_logz(const MatrixXl& logpot) {
    VectorXl r(logpot.cols());
    for (int i = 0; i < r.size(); ++i) r(i) = logsumexp(logpot.col(i));
    return r;
}

void MRFParameterizer::Pseudolikelihood::update_obj_score(lbfgsfloatval_t& fx, const MatrixXl& logpot, const VectorXl& logz, const size_t m, const string& seq, const float& sw) {
    for (auto it = param.v_offset.cbegin(); it != param.v_offset.cend(); ++it) {
        const size_t i = it->first;
        const int xi = data(m, i);
        if (xi < param.num_var) fx += -sw * (logpot(xi, i) - logz(i));
        else {
            float wi;
            const string s = param.abc.get_degeneracy(seq[i], &wi);
            for (auto c = s.cbegin(); c != s.cend(); ++c)
                fx += wi * (-sw * logpot(param.abc.get_idx(*c), i) + sw * logz(i));
        }
    }
}

void MRFParameterizer::Pseudolikelihood::update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const MatrixXl& logpot, const VectorXl& logz, const size_t m, const string& seq, const float& sw) {
    MatrixXl nodebel(param.num_var, param.num_node);
    for (size_t t = 0; t < param.num_var; ++t)
        nodebel.row(t) = (logpot.row(t) - logz.transpose()).unaryExpr(&exp);
    nodebel *= sw;
    MapMatrixXl dv(g, param.num_node, param.num_var);
    for (auto it = param.v_offset.cbegin(); it != param.v_offset.cend(); ++it) {
        const size_t i = it->first;
        int xi = data(m, i);
        if (xi < param.num_var) dv(i, xi) -= sw;
        else {
            float wi;
            const string s = param.abc.get_degeneracy(seq[i], &wi);
            for (auto c = s.cbegin(); c != s.cend(); ++c) dv(i, param.abc.get_idx(*c)) -= wi * sw;
        }
        dv.row(i) += nodebel.col(i);
    }
    for (auto it = param.w_offset.cbegin(); it != param.w_offset.cend(); ++it) {
        const int i = it->first.idx1;
        const int j = it->first.idx2;
        MapMatrixXl dw(g + it->second, param.num_var, param.num_var);
        int xi = data(m, i);
        int xj = data(m, j);
        if (xi < param.num_var && xj < param.num_var) {
            dw(xi, xj) -= sw * 2.;
            dw.col(xj) += nodebel.col(i);
            dw.row(xi) += nodebel.col(j);
        } else {
            float wi, wj, wt;
            const string si = param.abc.get_degeneracy(seq[i], &wi);
            const string sj = param.abc.get_degeneracy(seq[j], &wj);
            wt = wi * wj;
            for (auto ci = si.cbegin(); ci != si.cend(); ++ci) {
                xi = param.abc.get_idx(*ci);
                for (auto cj = sj.cbegin(); cj != sj.cend(); ++cj) {
                    xj = param.abc.get_idx(*cj);
                    dw(xi, xj) -= wt * sw * 2.;
                    for (int t = 0; t < param.num_var; ++t) {
                        dw.col(xj) += wt * nodebel.col(i);
                        dw.row(xi) += wt * nodebel.col(j);
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
  seq_weight(seq_weight),
  neff(neff),
  lp(NULL) {
    data = trace_to_var(traces, param.abc);
}

lbfgsfloatval_t MRFParameterizer::AsymPseudolikelihood::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        const float sw = seq_weight(i);
        const VectorXl logpot = calc_logpot(x, i, seq);
        const lbfgsfloatval_t logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, i, seq, sw);
        update_gradient(x, g, logpot, logz, i, seq, sw);
    }
    return fx;
}

VectorXl MRFParameterizer::AsymPseudolikelihood::calc_logpot(const lbfgsfloatval_t *x, const size_t m, const string& seq) {
    MapKVectorXl v(x, lp->num_var);
    VectorXl logpot(v);
    for (auto it = lp->w_offset.cbegin(); it != lp->w_offset.cend(); ++it) {
        const size_t i = it->first.idx1;
        const size_t j = it->first.idx2;
        MapKMatrixXl w(x + it->second, lp->num_var, lp->num_var);
        const int xi = data(m, i);
        const int xj = data(m, j);
        if (xj < lp->num_var) logpot += w.col(xj);
        else {
            float wj;
            const string sj = lp->abc.get_degeneracy(seq[j], &wj);
            for (auto c = sj.cbegin(); c != sj.cend(); ++c)
                logpot += wj * w.col(lp->abc.get_idx(*c));
        }
    }
    return logpot;
}

lbfgsfloatval_t MRFParameterizer::AsymPseudolikelihood::calc_logz(const VectorXl& logpot) {
    return logsumexp(logpot);
}

void MRFParameterizer::AsymPseudolikelihood::update_obj_score(lbfgsfloatval_t& fx, const VectorXl& logpot, const lbfgsfloatval_t& logz, const size_t m, const string& seq, const float& sw) {
    const size_t i = lp->get_obs_node();
    const int xi = data(m, i);
    if (xi < lp->num_var) fx += -sw * (logpot(xi) - logz);
    else {
        float wi;
        const string s = lp->abc.get_degeneracy(seq[i], &wi);
        for (auto c = s.cbegin(); c != s.cend(); ++c)
            fx += wi * (-sw * logpot(lp->abc.get_idx(*c)) + sw * logz);
    }
}

void MRFParameterizer::AsymPseudolikelihood::update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const VectorXl& logpot, const lbfgsfloatval_t& logz, const size_t m, const string& seq, const float& sw) {
    VectorXl nodebel(lp->num_var);
    nodebel = (logpot - VectorXl::Constant(lp->num_var, logz)).unaryExpr(&exp);
    nodebel *= sw;
    MapVectorXl dv(g + lp->v_beg_pos(), lp->num_var);
    const size_t i = lp->get_obs_node();
    int xi = data(m, i);
    if (xi < lp->num_var) dv(xi) -= sw;
    else {
        float wi;
        const string s = lp->abc.get_degeneracy(seq[i], &wi);
        for (auto c = s.cbegin(); c != s.cend(); ++c) dv(lp->abc.get_idx(*c)) -= wi * sw;
    }
    dv += nodebel;
    for (auto it = lp->w_offset.cbegin(); it != lp->w_offset.cend(); ++it) {
        const int i = it->first.idx1;
        const int j = it->first.idx2;
        MapMatrixXl dw(g + it->second, lp->num_var, lp->num_var);
        int xi = data(m, i);
        int xj = data(m, j);
        if (xi < lp->num_var && xj < lp->num_var) {
            dw(xi, xj) -= sw;
            dw.col(xj) += nodebel;
        } else {
            float wi, wj, wt;
            const string si = lp->abc.get_degeneracy(seq[i], &wi);
            const string sj = lp->abc.get_degeneracy(seq[j], &wj);
            wt = wi * wj;
            for (auto ci = si.cbegin(); ci != si.cend(); ++ci) {
                xi = lp->abc.get_idx(*ci);
                for (auto cj = sj.cbegin(); cj != sj.cend(); ++cj) {
                    xj = lp->abc.get_idx(*cj);
                    dw(xi, xj) -= wt * sw;
                    dw.col(xj) += wt * nodebel;
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
    for (auto it = funcs.cbegin(); it != funcs.cend(); ++it)
        fx += (*it)->evaluate(x, g, n, step);
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
    /* sequence weights */
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
            LocalParameter local_param(model, optim_opt, param, i);
            pll.set_local_param(&local_param);
            funcs.push_back(&pll);
            L2Regularization l2(local_param, opt.l2_opt);
            if (opt.regul == RegulMethod::REGUL_L2) funcs.push_back(&l2);
            ObjectiveFunction obj_func(funcs, local_param);
            int r = optimizer.optimize(&local_param, &obj_func);
            if (r < ret) ret = r;
            param.set_x(local_param);
            funcs.clear();
            //std::cout << "(obs_node, ret) = " << i << ", " << r << std::endl;
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
    SMMEmitProbEstimator emit_prob_estimator(ROBINSON_BGFREQ.get_array(aa), 
                                             BLOSUM62_MATRIX.get_array(aa), 
                                             3.2);
    PBSeqWeightEstimator seq_weight_estimator;
    ExpEntropyEffSeqNumEstimator eff_seq_num_estimator(aa, &seq_weight_estimator);
    TerminalGapRemover termi_gap_remover(aa, 0.1);
    ProfileBuilder profile_builder(aa);
    profile_builder.set_emit_prob_estimator(&emit_prob_estimator);
    profile_builder.set_seq_weight_estimator(&seq_weight_estimator);
    profile_builder.set_eff_seq_num_estimator(&eff_seq_num_estimator);
    profile_builder.set_msa_filter(&termi_gap_remover);
    Profile profile = profile_builder.build(traces);
    return profile.get_prob();
}

void MRFParameterizer::update_model(MRF& model, const SymmParameter& param) {
    for (size_t i = 0; i < model.get_length(); ++i)
        model.get_node(i).set_weight(MapKVectorXl(param.x + param.v_offset.at(i), param.num_var).cast<float>());
    const EdgeIndexVector eidxs = model.get_edge_idxs();
    for (auto it = eidxs.cbegin(); it != eidxs.cend(); ++it)
        model.get_edge(*it).set_weight(MapKMatrixXl(param.x + param.w_offset.at(*it), param.num_var, param.num_var).cast<float>());
}

void MRFParameterizer::update_model(MRF& model, const AsymParameter& param) {
    for (size_t i = 0; i < model.get_length(); ++i)
        model.get_node(i).set_weight(MapKVectorXl(param.x + param.v_offset.at(i), param.num_var).cast<float>());
    const EdgeIndexVector eidxs = model.get_edge_idxs();
    for (auto it = eidxs.cbegin(); it != eidxs.cend(); ++it) {
        MapKMatrixXl w1(param.x + param.w_offset.at(EdgeIndex(it->idx1, it->idx2)), param.num_var, param.num_var);
        MapKMatrixXl w2(param.x + param.w_offset.at(EdgeIndex(it->idx2, it->idx1)), param.num_var, param.num_var);
        model.get_edge(*it).set_weight((0.5 * (w1 + w2)).cast<float>());
    }
}

void MRFParameterizer::adjust_gauge(SymmParameter& param) {
    for (auto it = param.v_offset.cbegin(); it != param.v_offset.cend(); ++it) {
        MapVectorXl v(param.x + it->second, param.num_var);
        v = (v.array() - v.mean()).matrix();
    }
    for (auto it = param.w_offset.cbegin(); it != param.w_offset.cend(); ++it) {
        MapMatrixXl w(param.x + it->second, param.num_var, param.num_var);
        w = (w.array() - w.mean()).matrix();
    }
}

void MRFParameterizer::adjust_gauge(AsymParameter& param) {
    for (auto it = param.v_offset.cbegin(); it != param.v_offset.cend(); ++it) {
        MapVectorXl v(param.x + it->second, param.num_var);
        v = (v.array() - v.mean()).matrix();
        for (auto it2 = param.w_offset.cbegin(); it2 != param.w_offset.cend(); ++it2) {
            if (it2->first.idx1 != it->first) continue;
            MapMatrixXl w(param.x + it2->second, param.num_var, param.num_var);
            const VectorXl rowmn(w.rowwise().mean());
            const lbfgsfloatval_t totmn = w.mean();
            v += (rowmn.array() - totmn).matrix();
        }
    }
    for (auto it = param.w_offset.cbegin(); it != param.w_offset.cend(); ++it) {
        MapMatrixXl w(param.x + it->second, param.num_var, param.num_var);
        const lbfgsfloatval_t totmn = w.mean();
        const VectorXl rowmn(w.rowwise().mean());
        const VectorXl colmn(w.colwise().mean());
        for (size_t k = 0; k < param.num_var; ++k) {
            w.col(k) -= rowmn;
            w.row(k) -= colmn;
        }
        w = (w.array() + totmn).matrix();
    }
}
