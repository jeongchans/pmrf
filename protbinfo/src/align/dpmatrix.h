#ifndef _DPMATRIX_H_
#define _DPMATRIX_H_

#include <vector>

#include "common/numeric.h"
#include "dpbase.h"

using std::vector;

class DPMatrixVisitor;

class DPMatrixElement : public DPMatrixElementBase {
  public:
    virtual void accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col) = 0;
};

class MatchMatrixElement : public DPMatrixElement {
  public:
    virtual void accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class InsDelMatrixElement : public DPMatrixElement {
  public:
    virtual void accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class DelInsMatrixElement : public DPMatrixElement {
  public:
    virtual void accept(DPMatrixVisitor* visitor, const size_t& row, const size_t& col);
};

class DPMatrix {
  public:
    DPMatrix(const size_t& n, const size_t& m);

    vector<vector<MatchMatrixElement> > match;
    vector<vector<InsDelMatrixElement> > insdel;
    vector<vector<DelInsMatrixElement> > delins;
    const size_t rows, cols;

    void traverse(DPMatrixVisitor& visitor);
};

class DPMatrixVisitor {
  public:
    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col) = 0;
};

#endif
