#include "dpmatrix.h"

template <typename T>
void resize_2darray(vector<vector<T> >& arr, const size_t& m, const size_t& n) {
    arr.resize(m);
    for (size_t i = 0; i < m; ++i) arr[i].resize(n);
}

void MatchMatrixElement::accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_match(this, row, col);
}

void InsDelMatrixElement::accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_insdel(this, row, col);
}

void DelInsMatrixElement::accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_delins(this, row, col);
}

DPMatrix::DPMatrix(const size_t& n, const size_t& m) : rows(n + 1), cols(m + 1) {
    resize_2darray<MatchMatrixElement>(match, rows, cols);
    resize_2darray<InsDelMatrixElement>(insdel, rows, cols);
    resize_2darray<DelInsMatrixElement>(delins, rows, cols);
}

void DPMatrix::traverse(DPMatrixVisitor& visitor) {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            match[i][j].accept(&visitor, i, j);
            insdel[i][j].accept(&visitor, i, j);
            delins[i][j].accept(&visitor, i, j);
        }
    }
}
