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

double ExpEntropyEffSeqNumEstimator::estimate(const vector<string>& msa) const {
    int rows = msa.size();
    if (rows == 0) return 0.0;
    int cols = msa[0].size();
    VectorXf wt = seq_weight_estimator->estimate(msa);
    double s = 0;
    int n = abc.get_canonical_size(false);
    VectorXf p(n);
    for (int i = 0; i < cols; ++i) {
        p.setZero();
        for (int j = 0; j < rows; ++j)
            p += abc.get_count(msa[j][i]) * wt(j);
        p /= p.sum();
        for (int j = 0; j < n; ++j)
            if (p(j) > 0) s += p(j) * log(p(j));
    }
    return exp(-s / (double) cols);
}
