#include "parameterize.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "util/common.h"
#include "seq/profile.h"
#include "seq/bgfreq.h"
#include "seq/subsmat.h"

#include "core.h"

using blitz::pow2;
using blitz::sum;
using blitz::mean;
using blitz::pow;

/**
   @class MRFParameterizer::Parameter
 */

MRFParameterizer::Parameter::Parameter(const MRF& model, const Option& opt) : abc(model.get_alphabet()), length(model.get_length()) {
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    eidx.clear();
    int i = 0;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        eidx[EdgeIndex(pos->idx1, pos->idx2)] = i;
        i++;
    }
    num_var = model.get_num_var();
    n = length * num_var + eidx.size() * num_var * num_var;
    init();
    opt_param.linesearch = opt.linesearch;
    opt_param.delta = (lbfgsfloatval_t)(opt.delta);
    opt_param.past = opt.past;
    opt_param.max_iterations = opt.max_iterations;
}

int MRFParameterizer::Parameter::get_nidx(const int& i, const char& p) const {
    return i * num_var + abc.get_idx(p);
}

int MRFParameterizer::Parameter::get_eidx(const int& i, const int& j, const char& p, const char& q) const {
    int k = length * num_var;
    k += eidx.find(EdgeIndex(i, j))->second * num_var * num_var;
    k += abc.get_idx(p) * num_var;
    k += abc.get_idx(q);
    return k;
}

/**
   @class MRFParameterizer::L2Regularization
 */

void MRFParameterizer::L2Regularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    regularize_node(x, g, fx);
    regularize_edge(x, g, fx);
}

void MRFParameterizer::L2Regularization::regularize_node(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const double& lambda = opt.lambda1;
    string letters = param.abc.get_canonical();
    for (int i = 0; i < param.length; ++i) {
        for (int t = 0; t < param.num_var; ++t) {
            int k = param.get_nidx(i, letters[t]);
            fx += lambda * pow2(x[k]);
            g[k] += 2. * lambda * x[k];
        }
    }
}

void MRFParameterizer::L2Regularization::regularize_edge(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    double lambda = opt.lambda2;
    if (opt.sc) lambda *= 2. * ((double)(param.eidx.size())) / ((double)(param.length));
    string letters = param.abc.get_canonical();
    for (map<EdgeIndex, int>::const_iterator pos = param.eidx.begin(); pos != param.eidx.end(); ++pos) {
        int i = pos->first.idx1;
        int j = pos->first.idx2;
        for (int p = 0; p < param.num_var; ++p) {
            for (int q = 0; q < param.num_var; ++q) {
                int k = param.get_eidx(i, j, letters[p], letters[q]);
                fx += lambda * pow2(x[k]);
                g[k] += 2. * lambda * x[k];
            }
        }
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
    vector<string> msa = traces.get_matched_aseq_vec();
    msa = msa_analyzer.termi_gap_remover->filter(msa);
    seq_weight.resize(traces.size());
    seq_weight = msa_analyzer.seq_weight_estimator->estimate(msa);
    scale(seq_weight);
    seq_weight *= (double)msa_analyzer.eff_seq_num_estimator->estimate(msa);
}

lbfgsfloatval_t MRFParameterizer::ObjectiveFunction::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (int i = 0; i < param.n; ++i) g[i] = 0.;
    for (size_t i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_matched_aseq();
        double sw = seq_weight(i);
        Float2dArray logpot = calc_logpot(x, seq);
        Float1dArray logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, seq, sw);
        update_gradient(x, g, logpot, logz, seq, sw);
    }
    if (opt.regul == RegulMethod::RegulMethod::L2) l2_func.regularize(x, g, fx);
    return fx;
}

Float2dArray MRFParameterizer::ObjectiveFunction::calc_logpot(const lbfgsfloatval_t *x, const string& seq) {
    Float2dArray logpot = zeros(param.num_var, param.length);
    string letters = param.abc.get_canonical();
    FloatType wi, wj;
    for (int i = 0; i < param.length; ++i)
        for (int t = 0; t < param.num_var; ++t)
            logpot(t, i) += x[param.get_nidx(i, letters[t])];
    for (map<EdgeIndex, int>::const_iterator pos = param.eidx.begin(); pos != param.eidx.end(); ++pos) {
        int i = pos->first.idx1;
        int j = pos->first.idx2;
        string si = param.abc.get_degeneracy(seq[i], &wi);
        string sj = param.abc.get_degeneracy(seq[j], &wj);
        for (int t = 0; t < param.num_var; ++t) {
            for (string::iterator c = sj.begin(); c != sj.end(); ++c)
                logpot(t, i) += wj * x[param.get_eidx(i, j, letters[t], *c)];
            for (string::iterator c = si.begin(); c != si.end(); ++c)
                logpot(t, j) += wi * x[param.get_eidx(i, j, *c, letters[t])];
        }
    }
    return logpot;
}

Float1dArray MRFParameterizer::ObjectiveFunction::logsumexp(const Float2dArray& b) {
    Float1dArray B = row_max(b);
    Float1dArray r = zeros(b.rows());
    for (int j = 0; j < b.cols(); ++j)
        r += exp(b(ALL, j) - B);
    r = log(r) + B;
    return r;
}

Float1dArray MRFParameterizer::ObjectiveFunction::calc_logz(const Float2dArray& logpot) {
    Float1dArray logz = logsumexp(transpose(logpot));
    return logz;
}

void MRFParameterizer::ObjectiveFunction::update_obj_score(lbfgsfloatval_t& fx, const Float2dArray& logpot, const Float1dArray& logz, const string& seq, const double& sw) {
    for (int i = 0; i < param.length; ++i) {
        FloatType w;
        string s = param.abc.get_degeneracy(seq[i], &w);
        for (string::iterator c = s.begin(); c != s.end(); ++c)
            fx += w * (-sw * logpot(param.abc.get_idx(*c), i) + sw * logz(i));
    }
}

void MRFParameterizer::ObjectiveFunction::update_gradient(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const Float2dArray& logpot, const Float1dArray& logz, const string& seq, const double& sw) {
    string letters = param.abc.get_canonical();
    FloatType w, wi, wj;
    Float2dArray nodebel = zeros(param.num_var, param.length);
    for (int t = 0; t < param.num_var; ++t)
        nodebel(t, ALL) = exp(logpot(t, ALL) - logz);
    for (int i = 0; i < param.length; ++i) {
        string s = param.abc.get_degeneracy(seq[i], &w);
        for (string::iterator c = s.begin(); c != s.end(); ++c)
            g[param.get_nidx(i, *c)] -= w * sw;
        for (int t = 0; t < param.num_var; ++t)
            g[param.get_nidx(i, letters[t])] += sw * nodebel(t, i);
    }
    for (map<EdgeIndex, int>::const_iterator pos = param.eidx.begin(); pos != param.eidx.end(); ++pos) {
        int i = pos->first.idx1;
        int j = pos->first.idx2;
        string si = param.abc.get_degeneracy(seq[i], &wi);
        string sj = param.abc.get_degeneracy(seq[j], &wj);
        w = wi * wj;
        for (string::iterator ci = si.begin(); ci != si.end(); ++ci) {
            for (string::iterator cj = sj.begin(); cj != sj.end(); ++cj) {
                g[param.get_eidx(i, j, *ci, *cj)] -= w * sw * 2.;
                for (int t = 0; t < param.num_var; ++t) {
                    g[param.get_eidx(i, j, letters[t], *cj)] += w * sw * nodebel(t, i);
                    g[param.get_eidx(i, j, *ci, letters[t])] += w * sw * nodebel(t, j);
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
    Float2dArray psfm = zeros(param.length, param.num_var);
    FloatType gap_prob = 0.;
    psfm(ALL, Range(blitz::fromStart, 19)) = calc_profile(traces) * (1. - gap_prob);
    psfm(ALL, param.num_var - 1) = gap_prob;
    ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);
    LBFGS::Optimizer optimizer;
    int ret = optimizer.optimize(&param, &obj_func);
    model.set_psfm(psfm);
    update_model(model, param);
    return ret;
}

Float2dArray MRFParameterizer::calc_profile(const TraceVector& traces) {
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
        Float1dArray w = zeros(param.num_var);
        for (int t = 0; t < param.num_var; ++t) w(t) = (FloatType)param.x[param.get_nidx(i, letters[t])];
        model.get_node(i).set_weight(w);
    }
    for (map<EdgeIndex, int>::const_iterator pos = param.eidx.begin(); pos != param.eidx.end(); ++pos) {
        int i = pos->first.idx1;
        int j = pos->first.idx2;
        Float2dArray w = zeros(param.num_var, param.num_var);
        for (int p = 0; p < param.num_var; ++p)
            for (int q = 0; q < param.num_var; ++q) w(p, q) = (FloatType)param.x[param.get_eidx(i, j, letters[p], letters[q])];
        model.get_edge(i, j).set_weight(w);
    }
}
