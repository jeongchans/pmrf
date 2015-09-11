#include "hhdpmatrix.h"

template <typename T>
void resize_2darray(vector<vector<T> >& arr, const size_t& m, const size_t& n) {
    arr.resize(m);
    for (size_t i = 0; i < m; ++i) arr[i].resize(n);
}

void MMMatrixElement::accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_mm(this, row, col);
}

void MIMatrixElement::accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_mi(this, row, col);
}

void IMMatrixElement::accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_im(this, row, col);
}

void DGMatrixElement::accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_dg(this, row, col);
}

void GDMatrixElement::accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col) {
    visitor->visit_gd(this, row, col);
}

HHDPMatrix::HHDPMatrix(const size_t& n, const size_t& m) : rows(n + 1), cols(m + 1) {
    resize_2darray<MMMatrixElement>(mm, rows, cols);
    resize_2darray<MIMatrixElement>(mi, rows, cols);
    resize_2darray<IMMatrixElement>(im, rows, cols);
    resize_2darray<DGMatrixElement>(dg, rows, cols);
    resize_2darray<GDMatrixElement>(gd, rows, cols);
}

void HHDPMatrix::traverse(HHDPMatrixVisitor& visitor) {
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            mm[i][j].accept(&visitor, i, j);
            mi[i][j].accept(&visitor, i, j);
            im[i][j].accept(&visitor, i, j);
            dg[i][j].accept(&visitor, i, j);
            gd[i][j].accept(&visitor, i, j);
        }
    }
}
