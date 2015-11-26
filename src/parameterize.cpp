#include "parameterize.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

#include "util/common.h"
#include "seq/profile.h"
#include "core.h"
#include "seq/bgfreq.h"
#include "seq/subsmat.h"

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

///**
//   @class MRFParameterizer::NodePSSMRegularization
// */
//
//MRFParameterizer::NodePSSMRegularization::NodePSSMRegularization(const TraceVector& traces, Parameter& param, Option& opt) : param(param), opt(opt) {
//    if (opt.lambda != 0.) {
//        Float2dArray pssm = calc_pssm(traces);
//        string letters = param.abc.get_canonical();
//        mn.resize(param.length, param.num_var);
//        for (int i = 0; i < param.length; ++i) {
//            for (int t = 0; t < param.num_var; ++t) {
//                if (letters[t] == GAP_OPEN_SYMBOL) mn(i, t) = opt.gap_open;
//                else if (letters[t] == GAP_EXT_SYMBOL) mn(i, t) = opt.gap_ext;
//                else if (letters[t] == GAP_UNALI_SYMBOL) mn(i, t) = 0.;
//                else mn(i, t) = pssm(i, t);
//            }
//        }
//    }
//}
//
//Float2dArray MRFParameterizer::NodePSSMRegularization::calc_pssm(const TraceVector& traces) {
//    AminoAcid aa;
//    ProfileBuilder profile_builder(aa);
//    SMMEmitProbEstimator emit_prob_estimator(ROBINSON_BGFREQ.get_array(aa), 
//                                             BLOSUM62_MATRIX.get_array(aa), 
//                                             3.2);
//    profile_builder.set_emit_prob_estimator(&emit_prob_estimator);
//    PBSeqWeightEstimator seq_weight_estimator(aa);
//    profile_builder.set_seq_weight_estimator(&seq_weight_estimator);
//    ExpEntropyEffSeqNumEstimator eff_seq_num_estimator(aa, &seq_weight_estimator);
//    profile_builder.set_eff_seq_num_estimator(&eff_seq_num_estimator);
//    TerminalGapRemover termi_gap_remover(aa, 0.1);
//    profile_builder.set_msa_filter(&termi_gap_remover);
//    Profile profile = profile_builder.build(traces);
//    return profile.get_ll();
//}
//
//void MRFParameterizer::NodePSSMRegularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
//    const double& lambda = opt.lambda;
//    string letters = param.abc.get_canonical();
//    for (int i = 0; i < param.length; ++i) {
//        for (int t = 0; t < param.num_var; ++t) {
//            int k = param.get_nidx(i, letters[t]);
//            double d = x[k] - mn(i, t);
//            fx += lambda * pow2(d);
//            g[k] += 2. * lambda * (d);
//        }
//    }
//}

/**
   @class MRFParameterizer::NodeProfileRegularization
 */

MRFParameterizer::NodeProfileRegularization::NodeProfileRegularization(const TraceVector& traces, Parameter& param, Option& opt) : param(param), opt(opt) {
    if (opt.lambda != 0.) {
        Float2dArray freq = calc_profile(traces);
        string letters = param.abc.get_canonical();
        mn.resize(param.length, param.num_var);
        for (int i = 0; i < param.length; ++i) {
            for (int t = 0; t < param.num_var; ++t) {
                if (param.abc.is_gap(letters[t])) mn(i, t) = 0.;
                else mn(i, t) = log(freq(i, t) * (1. - opt.gap_prob) / opt.gap_prob);
            }
        }
    }
}

Float2dArray MRFParameterizer::NodeProfileRegularization::calc_profile(const TraceVector& traces) {
    AminoAcid aa;
    ProfileBuilder profile_builder(aa);
    SMMEmitProbEstimator emit_prob_estimator(ROBINSON_BGFREQ.get_array(aa), 
                                             BLOSUM62_MATRIX.get_array(aa), 
                                             3.2);
    profile_builder.set_emit_prob_estimator(&emit_prob_estimator);
    PBSeqWeightEstimator seq_weight_estimator(aa);
    profile_builder.set_seq_weight_estimator(&seq_weight_estimator);
    ExpEntropyEffSeqNumEstimator eff_seq_num_estimator(aa, &seq_weight_estimator);
    profile_builder.set_eff_seq_num_estimator(&eff_seq_num_estimator);
    TerminalGapRemover termi_gap_remover(aa, 0.1);
    profile_builder.set_msa_filter(&termi_gap_remover);
    Profile profile = profile_builder.build(traces);
    return profile.get_prob();
}

void MRFParameterizer::NodeProfileRegularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const double& lambda = opt.lambda;
    string letters = param.abc.get_canonical();
    for (int i = 0; i < param.length; ++i) {
        for (int t = 0; t < param.num_var; ++t) {
            int k = param.get_nidx(i, letters[t]);
            double d = x[k] - mn(i, t);
            fx += lambda * pow2(d);
            g[k] += 2. * lambda * (d);
        }
    }
}

/**
   @class MRFParameterizer::EdgeProfileRegularization
 */

MRFParameterizer::EdgeProfileRegularization::EdgeProfileRegularization(const TraceVector& traces, Parameter& param, Option& opt) : param(param), opt(opt) {
    if (opt.lambda != 0.) {
        if (opt.sc) lambda = 2. * ((double)(param.eidx.size())) / ((double)(param.length)) * opt.lambda;
        Float2dArray freq = calc_profile(traces);
        string letters = param.abc.get_canonical();
        mn.resize(param.length, param.num_var);
        for (int i = 0; i < param.length; ++i) {
            for (int t = 0; t < param.num_var; ++t) {
                if (param.abc.is_gap(letters[t])) mn(i, t) = 0.;
                else mn(i, t) = log(freq(i, t) * (1. - opt.gap_prob) / opt.gap_prob);
            }
        }
    }
}

Float2dArray MRFParameterizer::EdgeProfileRegularization::calc_profile(const TraceVector& traces) {
    AminoAcid aa;
    ProfileBuilder profile_builder(aa);
    SMMEmitProbEstimator emit_prob_estimator(ROBINSON_BGFREQ.get_array(aa), 
                                             BLOSUM62_MATRIX.get_array(aa), 
                                             3.2);
    profile_builder.set_emit_prob_estimator(&emit_prob_estimator);
    PBSeqWeightEstimator seq_weight_estimator(aa);
    profile_builder.set_seq_weight_estimator(&seq_weight_estimator);
    ExpEntropyEffSeqNumEstimator eff_seq_num_estimator(aa, &seq_weight_estimator);
    profile_builder.set_eff_seq_num_estimator(&eff_seq_num_estimator);
    TerminalGapRemover termi_gap_remover(aa, 0.1);
    profile_builder.set_msa_filter(&termi_gap_remover);
    Profile profile = profile_builder.build(traces);
    return profile.get_prob();
}

void MRFParameterizer::EdgeProfileRegularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    string letters = param.abc.get_canonical();
    for (map<EdgeIndex, int>::const_iterator pos = param.eidx.begin(); pos != param.eidx.end(); ++pos) {
        int i = pos->first.idx1;
        int j = pos->first.idx2;
        for (int p = 0; p < param.num_var; ++p) {
            for (int q = 0; q < param.num_var; ++q) {
                int k = param.get_eidx(i, j, letters[p], letters[q]);
                double d = x[k] - (mn(i, p) + mn(j, q));
                fx += lambda * pow2(d);
                g[k] += 2. * lambda * d;
            }
        }
    }
}

/**
   @class MRFParameterizer::NodeL2Regularization
 */

void MRFParameterizer::NodeL2Regularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
    const double& lambda = opt.lambda;
    string letters = param.abc.get_canonical();
    for (int i = 0; i < param.length; ++i) {
        for (int t = 0; t < param.num_var; ++t) {
            int k = param.get_nidx(i, letters[t]);
            fx += lambda * pow2(x[k]);
            g[k] += 2. * lambda * x[k];
        }
    }
}

/**
   @class MRFParameterizer::EdgeL2Regularization
 */

MRFParameterizer::EdgeL2Regularization::EdgeL2Regularization(Parameter& param, Option& opt) : param(param), opt(opt), lambda(opt.lambda) {
    if (opt.sc) lambda = 2. * ((double)(param.eidx.size())) / ((double)(param.length)) * opt.lambda;
}

void MRFParameterizer::EdgeL2Regularization::regularize(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, lbfgsfloatval_t& fx) {
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
  node_l2_func(param, opt.node_l2_opt),
  //node_pssm_func(traces, param, opt.node_pssm_opt),
  node_pb_func(traces, param, opt.node_pb_opt),
  edge_l2_func(param, opt.edge_l2_opt),
  edge_pb_func(traces, param, opt.edge_pb_opt),
  msa_analyzer(msa_analyzer) {
    vector<string> msa = traces.get_trimmed_aseq_vec();
    msa = msa_analyzer.termi_gap_remover->filter(msa);
    seq_weight.resize(traces.size());
    seq_weight = msa_analyzer.seq_weight_estimator->estimate(msa);
    scale(seq_weight);
    seq_weight *= (double)msa_analyzer.eff_seq_num_estimator->estimate(msa);
}

lbfgsfloatval_t MRFParameterizer::ObjectiveFunction::evaluate(const lbfgsfloatval_t *x, lbfgsfloatval_t *g, const int n, const lbfgsfloatval_t step) {
    lbfgsfloatval_t fx = 0.;
    for (int i = 0; i < param.n; ++i) g[i] = 0.;
    for (int i = 0; i < traces.size(); ++i) {
        string seq = traces[i].get_trimmed_aseq();
        double sw = seq_weight(i);
        Float2dArray logpot = calc_logpot(x, seq, sw);
        Float1dArray logz = calc_logz(logpot);
        update_obj_score(fx, logpot, logz, seq, sw);
        update_gradient(x, g, logpot, logz, seq, sw);
    }
    if (opt.node_regul == NodeRegulMethod::L2) node_l2_func.regularize(x, g, fx);
    //else if (opt.node_regul == NodeRegulMethod::PSSM) node_pssm_func.regularize(x, g, fx);
    else if (opt.node_regul == NodeRegulMethod::PROFILE) node_pb_func.regularize(x, g, fx);
    if (opt.edge_regul == EdgeRegulMethod::L2) edge_l2_func.regularize(x, g, fx);
    if (opt.edge_regul == EdgeRegulMethod::PROFILE) edge_pb_func.regularize(x, g, fx);
    return fx;
}

Float2dArray MRFParameterizer::ObjectiveFunction::calc_logpot(const lbfgsfloatval_t *x, const string& seq, const double& sw) {
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
    ObjectiveFunction obj_func(traces, param, opt, msa_analyzer);
    LBFGS::Optimizer optimizer;
    int ret = optimizer.optimize(&param, &obj_func);
    update_model(model, param);
    return ret;
}

void MRFParameterizer::update_model(MRF& model, Parameter& param) {
    string letters = param.abc.get_canonical();
    for (int i = 0; i < model.get_length(); ++i) {
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
