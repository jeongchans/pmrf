#include "hmmalign.h"

#include <cmath>

/**
 * Common functions
 */

inline double _calc_column_score(const HMMMatchState& match1, const HMMMatchState& match2, const Float1dArray& bgfreq) {
    Float1dArray x = match1.get_emit();
    Float1dArray y = match2.get_emit();
    if (sum(x) == 0) x = bgfreq;
    if (sum(y) == 0) y = bgfreq;
    return log2(sum(x * y / bgfreq));
}

inline double _calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to) {
    if (idx < 0 || idx >= (int)hmm.get_length()) return 0.;
    double prob;
    if (from == MATCH) prob = hmm.get_match(idx).get_transit_to(to);
    else if (from == INSERT) prob = hmm.get_insert(idx).get_transit_to(to);
    else if (from == DELETE) prob = hmm.get_delete(idx).get_transit_to(to);
    else return 0.;
    return log2(prob);
}

PairwiseAlignment _make_alignment(const ProfileHMM& hmm1, const ProfileHMM& hmm2, const DPTrace& tr) {
    const string x = hmm1.get_seq();
    const string y = hmm2.get_seq();
    PairwiseAlignment r;
    for (deque<pair<int, int> >::const_iterator pos = tr.res_trace.begin(); pos != tr.res_trace.end(); ++pos) {
        if (pos->first == GAP) r.alignment[0].push_back('-', GAP);
        else r.alignment[0].push_back(x[pos->first], pos->first);
        if (pos->second == GAP) r.alignment[1].push_back('-', GAP);
        else r.alignment[1].push_back(y[pos->second], pos->second);
    }
    return r;
}

/**
 * Local HMM-HMM alignment
 */

void LocalHHMatrixInitializer::visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) {
    if (row == 0 || col == 0) {
        element->score = 0;
        element->from = {NULL, -1, -1};
    } else {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    }
}

void LocalHHMatrixInitializer::visit_mi(MIMatrixElement* element, const size_t& row, const size_t&) {
    if (row == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else {
        element->score = 0;
        element->from = {NULL, -1, -1};
    }
}

void LocalHHMatrixInitializer::visit_im(IMMatrixElement* element, const size_t&, const size_t& col) {
    if (col == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else {
        element->score = 0;
        element->from = {NULL, -1, -1};
    }
}

void LocalHHMatrixInitializer::visit_dg(DGMatrixElement* element, const size_t& row, const size_t&) {
    if (row == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else {
        element->score = 0;
        element->from = {NULL, -1, -1};
    }
}

void LocalHHMatrixInitializer::visit_gd(GDMatrixElement* element, const size_t&, const size_t& col) {
    if (col == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else {
        element->score = 0;
        element->from = {NULL, -1, -1};
    }
}

void LocalHHMatrixOptimizer::visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0 && col > 0) {
        size_t x = row - 1;     // current residue index as well as previous matrix position
        size_t y = col - 1;
        int p = x - 1;  // previous residue index
        int q = y - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max5_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.mm[x][y]),
                make_pair(matrix.mi[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, INSERT, MATCH),
                          &matrix.mi[x][y]),
                make_pair(matrix.im[x][y].score +
                          calc_transit_score(hmm1, p, INSERT, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.im[x][y]),
                make_pair(matrix.dg[x][y].score +
                          calc_transit_score(hmm1, p, DELETE, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.dg[x][y]),
                make_pair(matrix.gd[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, DELETE, MATCH),
                          &matrix.gd[x][y]));
        double score = mx.first +
                       calc_column_score(hmm1.get_match(x), hmm2.get_match(y)) +
                       baseline;
        if (score >= 0) {
            element->score = score;
            element->from = {mx.second, (int)x, (int)y};
        } else {
            element->score = 0;
            element->from = {NULL, -1, -1};
        }
    }
}

void LocalHHMatrixOptimizer::visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        size_t x = row - 1;
        size_t y = col;
        int p = x - 1;  // previous residue index
        int q = y;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, INSERT),
                          &matrix.mm[x][y]),
                make_pair(matrix.mi[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, INSERT, INSERT),
                          &matrix.mi[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void LocalHHMatrixOptimizer::visit_im(IMMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        size_t x = row;
        size_t y = col - 1;
        int p = x;
        int q = y - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, INSERT) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.mm[x][y]),
                make_pair(matrix.im[x][y].score +
                          calc_transit_score(hmm1, p, INSERT, INSERT) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.im[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void LocalHHMatrixOptimizer::visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        size_t x = row - 1;
        size_t y = col;
        int p = x - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, DELETE),
                          &matrix.mm[x][y]),
                make_pair(matrix.dg[x][y].score +
                          calc_transit_score(hmm1, p, DELETE, DELETE),
                          &matrix.dg[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void LocalHHMatrixOptimizer::visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        size_t x = row;
        size_t y = col - 1;
        int q = y - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm2, q, MATCH, DELETE),
                          &matrix.mm[x][y]),
                make_pair(matrix.gd[x][y].score +
                          calc_transit_score(hmm2, q, DELETE, DELETE),
                          &matrix.gd[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

double LocalHHMatrixOptimizer::calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to) {
    return _calc_transit_score(hmm, idx, from, to);
}

double LocalHHMatrixOptimizer::calc_column_score(const HMMMatchState& match1, const HMMMatchState& match2) {
    return _calc_column_score(match1, match2, bgfreq);
}

pair<double, DPElementPosition> LocalHHMatrixOptimizer::find_optimum(const HHDPMatrix& matrix) {
    double opt_score = -INFINITY;
    DPElementPosition opt_pos = {NULL, -1, -1};
    pair<double, const HHDPMatrixElement*> p;
    for (size_t i = 0; i < matrix.rows; ++i) {
        for (size_t j = 0; j < matrix.cols; ++j) {
            p = make_pair(matrix.mm[i][j].score,  &matrix.mm[i][j]);
            if (opt_score < p.first) {
                opt_score = p.first;
                opt_pos = {p.second, (int)i, (int)j};
            }
        }
    }
    return make_pair(opt_score, opt_pos);
}

void LocalHHAlign::init_dp_align(const ProfileHMM& x, const ProfileHMM& y) {
    matrix = shared_ptr<HHDPMatrix>(
        new HHDPMatrix(x.get_length(), y.get_length()));
    initializer = shared_ptr<HHDPMatrixInitializer>(
        new LocalHHMatrixInitializer(x, y, *(matrix.get())));
    optimizer = shared_ptr<HHDPMatrixOptimizer>(
        new LocalHHMatrixOptimizer(bgfreq, x, y, baseline, *(matrix.get())));
}

PairwiseAlignment LocalHHAlign::make_alignment(const ProfileHMM& hmm1, const ProfileHMM& hmm2, const DPTrace& tr) {
    return _make_alignment(hmm1, hmm2, tr);
}

/**
 * Global HMM-HMM alignment
 */

void GlobalHHMatrixInitializer::visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) {
    if (row == 0 && col == 0) {
        element->score = 0;
        element->from = {NULL, -1, -1};
    } else {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    }
}

void GlobalHHMatrixInitializer::visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (row == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (row == 1) {
        element->score = 0;
        element->from = {&matrix.mm[row - 1][col], (int)row - 1, (int)col};
    } else {
        size_t x = row - 1;
        size_t y = col;
        int p = x - 1;
        element->score = matrix.mi[x][y].score + 
                         calc_transit_score(hmm1, p, MATCH, MATCH);
        element->from = {&matrix.mi[x][y], (int)x, (int)y};
    }
}

void GlobalHHMatrixInitializer::visit_im(IMMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (col == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (col == 1) {
        element->score = 0;
        element->from = {&matrix.mm[row][col - 1], (int)row, (int)col - 1};
    } else {
        size_t x = row;
        size_t y = col - 1;
        int q = y - 1;
        element->score = matrix.im[x][y].score + 
                         calc_transit_score(hmm2, q, MATCH, MATCH);
        element->from = {&matrix.im[x][y], (int)x, (int)y};
    }
}

void GlobalHHMatrixInitializer::visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (row == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (row == 1) {
        element->score = 0;
        element->from = {&matrix.mm[row - 1][col], (int)row - 1, (int)col};
    } else {
        size_t x = row - 1;
        size_t y = col;
        int p = x - 1;
        element->score = matrix.dg[x][y].score + 
                         calc_transit_score(hmm1, p, DELETE, DELETE);
        element->from = {&matrix.dg[x][y], (int)x, (int)y};
    }
}

void GlobalHHMatrixInitializer::visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (col == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (col == 1) {
        element->score = 0;
        element->from = {&matrix.mm[row][col - 1], (int)row, (int)col - 1};
    } else {
        size_t x = row;
        size_t y = col - 1;
        int q = y - 1;
        element->score = matrix.gd[x][y].score + 
                         calc_transit_score(hmm2, q, DELETE, DELETE);
        element->from = {&matrix.gd[x][y], (int)x, (int)y};
    }
}

double GlobalHHMatrixInitializer::calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to) {
    return _calc_transit_score(hmm, idx, from, to);
}

void GlobalHHMatrixOptimizer::visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0 && col > 0) {
        size_t x = row - 1;     // current residue index as well as previous matrix position
        size_t y = col - 1;
        int p = x - 1;  // previous residue index
        int q = y - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max5_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.mm[x][y]),
                make_pair(matrix.mi[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, INSERT, MATCH),
                          &matrix.mi[x][y]),
                make_pair(matrix.im[x][y].score +
                          calc_transit_score(hmm1, p, INSERT, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.im[x][y]),
                make_pair(matrix.dg[x][y].score +
                          calc_transit_score(hmm1, p, DELETE, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.dg[x][y]),
                make_pair(matrix.gd[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, DELETE, MATCH),
                          &matrix.gd[x][y]));
        double score = mx.first +
                       calc_column_score(hmm1.get_match(x), hmm2.get_match(y)) +
                       baseline;
        element->score = score;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void GlobalHHMatrixOptimizer::visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        size_t x = row - 1;
        size_t y = col;
        int p = x - 1;  // previous residue index
        int q = y;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, MATCH, INSERT),
                          &matrix.mm[x][y]),
                make_pair(matrix.mi[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, MATCH) +
                          calc_transit_score(hmm2, q, INSERT, INSERT),
                          &matrix.mi[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void GlobalHHMatrixOptimizer::visit_im(IMMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        size_t x = row;
        size_t y = col - 1;
        int p = x;
        int q = y - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, INSERT) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.mm[x][y]),
                make_pair(matrix.im[x][y].score +
                          calc_transit_score(hmm1, p, INSERT, INSERT) +
                          calc_transit_score(hmm2, q, MATCH, MATCH),
                          &matrix.im[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void GlobalHHMatrixOptimizer::visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        size_t x = row - 1;
        size_t y = col;
        int p = x - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm1, p, MATCH, DELETE),
                          &matrix.mm[x][y]),
                make_pair(matrix.dg[x][y].score +
                          calc_transit_score(hmm1, p, DELETE, DELETE),
                          &matrix.dg[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

void GlobalHHMatrixOptimizer::visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        size_t x = row;
        size_t y = col - 1;
        int q = y - 1;
        pair<double, const HHDPMatrixElement*> mx =
            max_pair<double, const HHDPMatrixElement*>(
                make_pair(matrix.mm[x][y].score +
                          calc_transit_score(hmm2, q, MATCH, DELETE),
                          &matrix.mm[x][y]),
                make_pair(matrix.gd[x][y].score +
                          calc_transit_score(hmm2, q, DELETE, DELETE),
                          &matrix.gd[x][y]));
        element->score = mx.first;
        element->from = {mx.second, (int)x, (int)y};
    }
}

double GlobalHHMatrixOptimizer::calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to) {
    return _calc_transit_score(hmm, idx, from, to);
}

double GlobalHHMatrixOptimizer::calc_column_score(const HMMMatchState& match1, const HMMMatchState& match2) {
    return _calc_column_score(match1, match2, bgfreq);
}

pair<double, DPElementPosition> GlobalHHMatrixOptimizer::find_optimum(const HHDPMatrix& matrix) {
    size_t i = matrix.rows - 1;
    size_t j = matrix.cols - 1;
    pair<double, const HHDPMatrixElement*> p =
        max5_pair<double, const HHDPMatrixElement*>(
            make_pair(matrix.mm[i][j].score,  &matrix.mm[i][j]),
            make_pair(matrix.mi[i][j].score,  &matrix.mi[i][j]),
            make_pair(matrix.im[i][j].score,  &matrix.im[i][j]),
            make_pair(matrix.dg[i][j].score,  &matrix.dg[i][j]),
            make_pair(matrix.gd[i][j].score,  &matrix.gd[i][j]));
    double opt_score = p.first;
    DPElementPosition opt_pos = {p.second, (int)i, (int)j};
    return make_pair(opt_score, opt_pos);
}

void GlobalHHAlign::init_dp_align(const ProfileHMM& x, const ProfileHMM& y) {
    matrix = shared_ptr<HHDPMatrix>(
        new HHDPMatrix(x.get_length(), y.get_length()));
    initializer = shared_ptr<HHDPMatrixInitializer>(
        new GlobalHHMatrixInitializer(x, y, *(matrix.get())));
    optimizer = shared_ptr<HHDPMatrixOptimizer>(
        new GlobalHHMatrixOptimizer(bgfreq, x, y, baseline, *(matrix.get())));
}

PairwiseAlignment GlobalHHAlign::make_alignment(const ProfileHMM& hmm1, const ProfileHMM& hmm2, const DPTrace& tr) {
    return _make_alignment(hmm1, hmm2, tr);
}
