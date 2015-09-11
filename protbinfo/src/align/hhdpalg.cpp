#include "hhdpalg.h"

#include <cmath>

using std::make_pair;

pair<double, DPElementPosition> HHDPMatrixOptimizer::optimize(HHDPMatrix& matrix) {
    matrix.traverse(*this);
    return find_optimum(matrix);
}

PairwiseAlignment HHDPAlgorithm::operator()(const ProfileHMM& x, const ProfileHMM& y) {
    init_dp_align(x, y);
    initializer->initialize(*matrix);
    pair<double, DPElementPosition> p = optimizer->optimize(*matrix);
    DPTrace tr = traceback(p.second);
    PairwiseAlignment r = make_alignment(x, y, tr);
    r.property["score"] = p.first;
    return r;
}

DPTrace HHDPAlgorithm::traceback(const DPElementPosition& beg) {
    DPTrace tr;
    const DPElementPosition *curr = &beg;
    while (curr->element != NULL) {
        const DPElementPosition *prev = &(curr->element->from);
        if (prev->element == NULL) break;
        int row = curr->row;
        int col = curr->col;
        if (&(matrix->mm[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(row - 1, col - 1));
        } else if (&(matrix->mi[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(row - 1, GAP));
        } else if (&(matrix->im[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(GAP, col - 1));
        } else if (&(matrix->dg[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(row - 1, GAP));
        } else if (&(matrix->gd[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(GAP, col - 1));
        }
        tr.element_trace.push_front(*curr);
        curr = prev;
    }
    return tr;
}
