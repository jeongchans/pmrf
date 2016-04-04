#include "seqweight.h"

#include <cctype>

VectorXf PBSeqWeightEstimator::estimate(const vector<string>& msa) const {
    int rows = msa.size();
    int cols = msa[0].size();
    VectorXf seq_wt = VectorXf::Zero(rows); // sequence weights
    VectorXf alen = VectorXf::Zero(rows);   // number of residues for normalization
    char c;
    for (int j = 0; j < cols; ++j) {
        VectorXf res_wt(dim);  // residue weight
        res_wt = calc_residue_weight(msa, j);
        for (int i = 0; i < rows; ++i) {
            c = msa[i][j];
            if (is_allowed(c)) {
                seq_wt(i) += res_wt(abc_idx(c));
                alen(i) += 1;
            }
        }
    }
    for (int i = 0; i < rows; ++i)
        if (alen(i) > 0) seq_wt(i) /= alen(i);     // normalized by length
    if (seq_wt.sum() == 0) seq_wt.setOnes();
    return seq_wt / seq_wt.sum();
}

VectorXf PBSeqWeightEstimator::calc_residue_weight(const vector<string>& msa, const int& idx) const {
    VectorXf s = VectorXf::Zero(dim);  // number of times a particular residue appears
    char c;
    for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos) {
        c = (*pos)[idx];
        if (is_allowed(c)) s(abc_idx(c)) += 1;
    }
    double r = (double) (s.array() > 0).count();    // number of different residues
    VectorXf wt(dim);
    for (int k = 0; k < dim; ++k) {
        if (s(k) > 0) wt(k) = 1. / (r * s(k));
        else wt(k) = 0;
    }
    return wt;
}

bool PBSeqWeightEstimator::is_allowed(const char& c) const {
    if (isalpha(c)) return true;
    else return false;
}

int PBSeqWeightEstimator::abc_idx(const char& c) const {
    return toupper((int) c) - 'A';
}

VectorXf PSIBLASTPBSeqWeightEstimator::estimate(const vector<string>& msa) const {
    int rows = msa.size();
    int cols = msa[0].size();
    VectorXf seq_wt = VectorXf::Zero(rows); // sequence weights
    char c;
    for (int j = 0; j < cols; ++j) {
        VectorXf res_wt(dim);  // residue weight
        res_wt = calc_residue_weight(msa, j);
        if (res_wt.sum() == 0) continue;
        for (int i = 0; i < rows; ++i) {
            c = msa[i][j];
            seq_wt(i) += res_wt(abc_idx(c));
        }
    }
    if (seq_wt.sum() == 0) seq_wt.setOnes();
    return seq_wt / seq_wt.sum();
}

VectorXf PSIBLASTPBSeqWeightEstimator::calc_residue_weight(const vector<string>& msa, const int& idx) const {
    VectorXf s = VectorXf::Zero(dim);  // number of times a particular residue appears
    char c;
    for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos) {
        c = (*pos)[idx];
        s(abc_idx(c)) += 1;
    }
    double r = (double) (s.array() > 0).count();    // number of different residues
    VectorXf wt(dim);
    if (r == 1) {
        wt.setZero();
        return wt;
    }
    for (int k = 0; k < dim; ++k) {
        if (s(k) > 0) wt(k) = 1. / (r * s(k));
        else wt(k) = 0;
    }
    return wt;
}

int PSIBLASTPBSeqWeightEstimator::abc_idx(const char& c) const {
    if (isalpha(c)) return toupper((int) c) - 'A';
    else return gap_idx;
}
