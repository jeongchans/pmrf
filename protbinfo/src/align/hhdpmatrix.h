#ifndef _HHDPMATRIX_H_
#define _HHDPMATRIX_H_

#include <vector>

#include "common/numeric.h"
#include "dpbase.h"

using std::vector;

class HHDPMatrixVisitor;

class HHDPMatrixElement : public DPMatrixElementBase {
  public:
    virtual void accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col) = 0;
};

class MMMatrixElement : public HHDPMatrixElement {
  public:
    virtual void accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class MIMatrixElement : public HHDPMatrixElement {
  public:
    virtual void accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class IMMatrixElement : public HHDPMatrixElement {
  public:
    virtual void accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class DGMatrixElement : public HHDPMatrixElement {
  public:
    virtual void accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class GDMatrixElement : public HHDPMatrixElement {
  public:
    virtual void accept(HHDPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class HHDPMatrix {
  public:
    HHDPMatrix(const size_t& n, const size_t& m);

    vector<vector<MMMatrixElement> > mm;
    vector<vector<MIMatrixElement> > mi;
    vector<vector<IMMatrixElement> > im;
    vector<vector<DGMatrixElement> > dg;
    vector<vector<GDMatrixElement> > gd;
    const size_t rows, cols;

    void traverse(HHDPMatrixVisitor& visitor);
};

class HHDPMatrixVisitor {
  public:
    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col) = 0;
};

#endif

