#include "dpalg.h"

#include <cmath>

using std::make_pair;

pair<double, DPElementPosition> DPMatrixOptimizer::optimize(DPMatrix& matrix) {
    matrix.traverse(*this);
    return find_optimum(matrix);
}

PairwiseAlignment DPAlgorithm::operator()(const string& x, const string& y) {
    init_dp_align(x, y);
    initializer->initialize(*matrix);
    pair<double, DPElementPosition> p = optimizer->optimize(*matrix);
    DPTrace tr = traceback(p.second);
    PairwiseAlignment r = make_alignment(x, y, tr);
    r.property["score"] = p.first;
    return r;
}

DPTrace DPAlgorithm::traceback(const DPElementPosition& beg) {
    DPTrace tr;
    const DPElementPosition *curr = &beg;
    while (curr->element != NULL) {
        const DPElementPosition *prev = &(curr->element->from);
        if (prev->element == NULL) break;
        int row = curr->row;
        int col = curr->col;
        if (&(matrix->match[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(row - 1, col - 1));
        } else if (&(matrix->insdel[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(row - 1, GAP));
        } else if (&(matrix->delins[row][col]) == curr->element) {
            tr.res_trace.push_front(make_pair(GAP, col - 1));
        }
        tr.element_trace.push_front(*curr);
        curr = prev;
    }
    return tr;
}
