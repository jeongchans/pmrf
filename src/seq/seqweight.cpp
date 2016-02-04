#include "seqweight.h"

#include <cctype>

Float1dArray PBSeqWeightEstimator::estimate(const vector<string>& msa) const {
    int rows = msa.size();
    int cols = msa[0].size();
    Float1dArray seq_wt(rows); // sequence weights
    seq_wt = 0;
    Float1dArray alen(rows);   // number of residues for normalization
    alen = 0;
    char c;
    for (int j = 0; j < cols; ++j) {
        Float1dArray res_wt(dim);  // residue weight
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
    if (blitz::sum(seq_wt) == 0) seq_wt = 1; 
    scale(seq_wt);
    return seq_wt;
}

Float1dArray PBSeqWeightEstimator::calc_residue_weight(const vector<string>& msa, const int& idx) const {
    Float1dArray s(dim);  // number of times a particular residue appears
    s = 0;
    char c;
    for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos) {
        c = (*pos)[idx];
        if (is_allowed(c)) s(abc_idx(c)) += 1;
    }
    double r = (double) blitz::count(s > 0);    // number of different residues
    Float1dArray wt(dim);
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

Float1dArray PSIBLASTPBSeqWeightEstimator::estimate(const vector<string>& msa) const {
    int rows = msa.size();
    int cols = msa[0].size();
    Float1dArray seq_wt(rows); // sequence weights
    seq_wt = 0;
    char c;
    for (int j = 0; j < cols; ++j) {
        Float1dArray res_wt(dim);  // residue weight
        res_wt = calc_residue_weight(msa, j);
        if (blitz::sum(res_wt) == 0) continue;
        for (int i = 0; i < rows; ++i) {
            c = msa[i][j];
            seq_wt(i) += res_wt(abc_idx(c));
        }
    }
    if (blitz::sum(seq_wt) == 0) seq_wt = 1; 
    scale(seq_wt);
    return seq_wt;
}

Float1dArray PSIBLASTPBSeqWeightEstimator::calc_residue_weight(const vector<string>& msa, const int& idx) const {
    Float1dArray s(dim);  // number of times a particular residue appears
    s = 0;
    char c;
    for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos) {
        c = (*pos)[idx];
        s(abc_idx(c)) += 1;
    }
    double r = (double) blitz::count(s > 0);    // number of different residues
    Float1dArray wt(dim);
    if (r == 1) {
        wt = 0;
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