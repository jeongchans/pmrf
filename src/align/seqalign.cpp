#include "seqalign.h"

#include <cmath>

/**
 * Needleman-Wunch
 */

void NWMatrixInitializer::visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) {
    if (row == 0 && col == 0) {
        element->score = 0;
        element->from = {NULL, -1, -1};
    } else {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    }
}

void NWMatrixInitializer::visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col) {
    if (row == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (row == 1) {
        element->score = - gap_open;
        element->from = {&matrix.match[row - 1][col], (int)row - 1, (int)col};
    } else {
        element->score = - gap_open - (row - 1) * gap_ext;
        element->from = {&matrix.insdel[row - 1][col], (int)row - 1, (int)col};
    }
}

void NWMatrixInitializer::visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col) {
    if (col == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else if (col == 1) {
        element->score = - gap_open;
        element->from = {&matrix.match[row][col - 1], (int)row, (int)col - 1};
    } else {
        element->score = - gap_open - (col - 1) * gap_ext;
        element->from = {&matrix.delins[row][col - 1], (int)row, (int)col - 1};
    }
}

void NWMatrixOptimizer::visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0 && col > 0) {
        pair<double, const DPMatrixElement*> p = 
            max3_pair<double, const DPMatrixElement*>(
                make_pair(matrix.match[row - 1][col - 1].score, 
                          &matrix.match[row - 1][col - 1]),
                make_pair(matrix.insdel[row - 1][col - 1].score, 
                          &matrix.insdel[row - 1][col - 1]),
                make_pair(matrix.delins[row - 1][col - 1].score, 
                          &matrix.delins[row - 1][col - 1]));
        element->score = p.first + 
                         subsmat.get_value(seq1[row - 1], seq2[col - 1]) +
                         baseline;
        element->from = {p.second, (int)row - 1, (int)col - 1};
    }
}

void NWMatrixOptimizer::visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        pair<double, const DPMatrixElement*> p = 
            max_pair<double, const DPMatrixElement*>(
                make_pair(matrix.match[row - 1][col].score - gap_open,
                          &matrix.match[row - 1][col]),
                make_pair(matrix.insdel[row - 1][col].score - gap_ext,
                          &matrix.insdel[row - 1][col]));
        element->score = p.first;
        element->from = {p.second, (int)row - 1, (int)col};
    }
}

void NWMatrixOptimizer::visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        pair<double, const DPMatrixElement*> p = 
            max_pair<double, const DPMatrixElement*>(
                make_pair(matrix.match[row][col - 1].score - gap_open,
                          &matrix.match[row][col - 1]),
                make_pair(matrix.delins[row][col - 1].score - gap_ext,
                          &matrix.delins[row][col - 1]));
        element->score = p.first;
        element->from = {p.second, (int)row, (int)col - 1};
    }
}

pair<double, DPElementPosition> NWMatrixOptimizer::find_optimum(const DPMatrix& matrix) {
    size_t i = matrix.rows - 1;
    size_t j = matrix.cols - 1;
    pair<double, const DPMatrixElement*> p =
        max3_pair<double, const DPMatrixElement*>(
            make_pair(matrix.match[i][j].score,  &matrix.match[i][j]),
            make_pair(matrix.insdel[i][j].score, &matrix.insdel[i][j]),
            make_pair(matrix.delins[i][j].score, &matrix.delins[i][j]));
    double opt_score = p.first;
    DPElementPosition opt_pos = {p.second, (int)i, (int)j};
    return make_pair(opt_score, opt_pos);
}

void NWAlign::init_dp_align(const string& x, const string& y) {
    matrix = shared_ptr<DPMatrix>(
        new DPMatrix(x.size(), y.size()));
    initializer = shared_ptr<DPMatrixInitializer>(
        new NWMatrixInitializer(gap_open, gap_ext, *(matrix.get())));
    optimizer = shared_ptr<DPMatrixOptimizer>(
        new NWMatrixOptimizer(subsmat, x, y, gap_open, gap_ext, baseline, *(matrix.get())));
}

PairwiseAlignment NWAlign::make_alignment(const string& x, const string& y, const DPTrace& tr) {
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
 * Smith-Waterman
 */

void SWMatrixInitializer::visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) {
    if (row == 0 && col == 0) {
        element->score = 0;
        element->from = {NULL, -1, -1};
    } else {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    }
}

void SWMatrixInitializer::visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t&) {
    if (row == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else {
        element->score = 0;
        element->from = {NULL, -1, -1};
    }
}

void SWMatrixInitializer::visit_delins(DelInsMatrixElement* element, const size_t&, const size_t& col) {
    if (col == 0) {
        element->score = -INFINITY;
        element->from = {NULL, -1, -1};
    } else {
        element->score = 0;
        element->from = {NULL, -1, -1};
    }
}

void SWMatrixOptimizer::visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0 && col > 0) {
        pair<double, const DPMatrixElement*> p = 
            max3_pair<double, const DPMatrixElement*>(
                make_pair(matrix.match[row - 1][col - 1].score, 
                          &matrix.match[row - 1][col - 1]),
                make_pair(matrix.insdel[row - 1][col - 1].score, 
                          &matrix.insdel[row - 1][col - 1]),
                make_pair(matrix.delins[row - 1][col - 1].score, 
                          &matrix.delins[row - 1][col - 1]));
        double score = p.first + 
                       subsmat.get_value(seq1[row - 1], seq2[col - 1]) +
                       baseline;
        if (score >= 0) {
            element->score = score;
            element->from = {p.second, (int)row - 1, (int)col - 1};
        } else {
            element->score = 0;
            element->from = {NULL, -1, -1};
        }
    }
}

void SWMatrixOptimizer::visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col) {
    if (row > 0) {
        pair<double, const DPMatrixElement*> p = 
            max_pair<double, const DPMatrixElement*>(
                make_pair(matrix.match[row - 1][col].score - gap_open,
                          &matrix.match[row - 1][col]),
                make_pair(matrix.insdel[row - 1][col].score - gap_ext,
                          &matrix.insdel[row - 1][col]));
        element->score = p.first;
        element->from = {p.second, (int)row - 1, (int)col};
    }
}

void SWMatrixOptimizer::visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col) {
    if (col > 0) {
        pair<double, const DPMatrixElement*> p = 
            max_pair<double, const DPMatrixElement*>(
                make_pair(matrix.match[row][col - 1].score - gap_open,
                          &matrix.match[row][col - 1]),
                make_pair(matrix.delins[row][col - 1].score - gap_ext,
                          &matrix.delins[row][col - 1]));
        element->score = p.first;
        element->from = {p.second, (int)row, (int)col - 1};
    }
}

pair<double, DPElementPosition> SWMatrixOptimizer::find_optimum(const DPMatrix& matrix) {
    double opt_score = -INFINITY;
    DPElementPosition opt_pos = {NULL, -1, -1};
    pair<double, const DPMatrixElement*> p;
    for (size_t i = 0; i < matrix.rows; ++i) {
        for (size_t j = 0; j < matrix.cols; ++j) {
            p = max3_pair<double, const DPMatrixElement*>(
                    make_pair(matrix.match[i][j].score,  &matrix.match[i][j]),
                    make_pair(matrix.insdel[i][j].score, &matrix.insdel[i][j]),
                    make_pair(matrix.delins[i][j].score, &matrix.delins[i][j]));
            if (opt_score < p.first) {
                opt_score = p.first;
                opt_pos = {p.second, (int)i, (int)j};
            }
        }
    }
    return make_pair(opt_score, opt_pos);
}

void SWAlign::init_dp_align(const string& x, const string& y) {
    matrix = shared_ptr<DPMatrix>(
        new DPMatrix(x.size(), y.size()));
    initializer = shared_ptr<DPMatrixInitializer>(
        new SWMatrixInitializer(gap_open, gap_ext, *(matrix.get())));
    optimizer = shared_ptr<DPMatrixOptimizer>(
        new SWMatrixOptimizer(subsmat, x, y, gap_open, gap_ext, baseline, *(matrix.get())));
}

PairwiseAlignment SWAlign::make_alignment(const string& x, const string& y, const DPTrace& tr) {
    PairwiseAlignment r;
    for (deque<pair<int, int> >::const_iterator pos = tr.res_trace.begin(); pos != tr.res_trace.end(); ++pos) {
        if (pos->first == GAP) r.alignment[0].push_back('-', GAP);
        else r.alignment[0].push_back(x[pos->first], pos->first);
        if (pos->second == GAP) r.alignment[1].push_back('-', GAP);
        else r.alignment[1].push_back(y[pos->second], pos->second);
    }
    return r;
}
