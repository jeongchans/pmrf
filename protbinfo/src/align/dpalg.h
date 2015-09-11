#ifndef _DPALG_H_
#define _DPALG_H_

#include <memory>
#include <string>

#include "dpmatrix.h"
#include "alignment.h"

using std::pair;
using std::shared_ptr;
using std::string;

class DPMatrixInitializer : public DPMatrixVisitor {
  public:
    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col) = 0;

    void initialize(DPMatrix& matrix) { matrix.traverse(*this); }
};

class DPMatrixOptimizer : public DPMatrixVisitor {
  public:
    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col) = 0;

    pair<double, DPElementPosition> optimize(DPMatrix& matrix);

  protected:
    virtual pair<double, DPElementPosition> find_optimum(const DPMatrix& matrix) = 0;
};

class DPAlgorithm {
  public:
    PairwiseAlignment operator()(const string& x, const string& y);

  protected:
    shared_ptr<DPMatrix> matrix;
    shared_ptr<DPMatrixInitializer> initializer;
    shared_ptr<DPMatrixOptimizer> optimizer;

    DPTrace traceback(const DPElementPosition& beg);

    virtual void init_dp_align(const string& x, const string& y) = 0;
    virtual PairwiseAlignment make_alignment(const string& x, const string& y, const DPTrace& tr) = 0;
};

#endif
