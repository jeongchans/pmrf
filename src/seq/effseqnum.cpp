#include "effseqnum.h"

double RTEffSeqNumEstimator::estimate(const vector<string>& msa) const {
    int cols = msa[0].size();
    VectorXf num_res_types(cols);  // number of distinct residue types at each position
    char c;
    for (int i = 0; i < cols; ++i) {
        string column = "";
        for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos) {
            c = (*pos)[i];
            if (abc.is_canonical(c, false)) column += c;
        }
        num_res_types(i) = calc_num_res_type(column);
    }
    return num_res_types.mean();
}

size_t RTEffSeqNumEstimator::calc_num_res_type(const string& column) const {
    VectorXf n = VectorXf::Zero(abc.get_canonical_size(false));
    for (string::const_iterator pos = column.begin(); pos != column.end(); ++pos)
        n(abc.get_idx(*pos)) += 1;
    return (n.array() != 0).count();
}

double ExpEntropyEffSeqNumEstimator::estimate(const vector<string>& msa, const VectorXf& wt) const {
    int rows = msa.size();
    if (rows == 0) return 0.;
    int cols = msa[0].size();
    double s = 0;
    int m = 0;
    int n = abc.get_canonical_size(false);
    for (int i = 0; i < cols; ++i) {
        VectorXf p = VectorXf::Zero(n);
        for (int j = 0; j < rows; ++j)
            p += wt(j) * abc.get_count(msa[j][i]);
        if (p.sum() > 0) {
            p /= p.sum();
            s += -(p.array() > 0).select(p.cwiseProduct(p.unaryExpr(&log)), 0.).sum();
            ++m;
        }
    }
    if (m > 0) return exp(s / (double) m);
    else return 0.;
}

double ExpJointEntropyEffSeqNumEstimator::estimate(const vector<string>& msa, const VectorXf& wt) const {
    int rows = msa.size();
    if (rows == 0) return 0.0;
    int cols = msa[0].size();
    double s = 0;
    int n = abc.get_canonical_size(false);
    for (int i = 0; i < cols; ++i) {
        for (int j = i + 1; j < cols; ++j) {
            MatrixXf p = MatrixXf::Zero(n, n);
            for (int m = 0; m < rows; ++m)
                p += wt(m) * abc.get_count(msa[m][i]) * abc.get_count(msa[m][j]).transpose();
            p /= p.sum();
            s += -(p.array() > 0).select(p.cwiseProduct(p.unaryExpr(&log)), 0.).sum();
        }
    }
    return exp(s / (double) (cols * (cols - 1) / 2));
}

double ClstrEffSeqNumEstimator::estimate(const vector<string>& msa) const {
    ClstrSeqWeightEstimator clstr_sw_estimator(abc, maxidt);
    return clstr_sw_estimator.estimate(msa).sum();
}
