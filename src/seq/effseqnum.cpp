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
    if (rows == 0) return 0.0;
    int cols = msa[0].size();
    double s = 0;
    int n = abc.get_canonical_size(false);
    for (int i = 0; i < cols; ++i) {
        VectorXf p = VectorXf::Zero(n);
        for (int j = 0; j < rows; ++j)
            p += wt(j) * abc.get_count(msa[j][i]);
        p /= p.sum();
        s += -(p.array() > 0).select(p.cwiseProduct(p.unaryExpr(&log)), 0.).sum();
    }
    return exp(s / (double) cols);
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
    size_t rows = msa.size();
    VectorXf memnum = VectorXf::Ones(rows);
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = i + 1; j < rows; ++j) {
            if (calc_identity(msa[i], msa[j]) > maxidt) {
                memnum(i) += 1;
                memnum(j) += 1;
            }
        }
    }
    return VectorXf::Ones(rows).cwiseQuotient(memnum).sum();
}

float ClstrEffSeqNumEstimator::calc_identity(const string& seq1, const string& seq2) const {
    size_t n = seq1.size();
    float f = 0.;
    for (size_t i = 0; i < n; ++i)
        if (seq1[i] == seq2[i]) ++f;
    return f / (float) n;
}
