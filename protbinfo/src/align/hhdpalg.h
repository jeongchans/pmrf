#ifndef _HHDPALG_H_
#define _HHDPALG_H_

#include <memory>
#include <string>

#include "hhdpmatrix.h"
#include "alignment.h"
#include "msaanalysis/profilehmm.h"

using std::pair;
using std::shared_ptr;
using std::string;

class HHDPMatrixInitializer : public HHDPMatrixVisitor {
  public:
    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col) = 0;

    void initialize(HHDPMatrix& matrix) { matrix.traverse(*this); }
};

class HHDPMatrixOptimizer : public HHDPMatrixVisitor {
  public:
    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col) = 0;
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col) = 0;

    pair<double, DPElementPosition> optimize(HHDPMatrix& matrix);

  protected:
    virtual pair<double, DPElementPosition> find_optimum(const HHDPMatrix& matrix) = 0;
};

class HHDPAlgorithm {
  public:
    PairwiseAlignment operator()(const ProfileHMM& x, const ProfileHMM& y);

  protected:
    shared_ptr<HHDPMatrix> matrix;
    shared_ptr<HHDPMatrixInitializer> initializer;
    shared_ptr<HHDPMatrixOptimizer> optimizer;

    DPTrace traceback(const DPElementPosition& beg);

    virtual void init_dp_align(const ProfileHMM& x, const ProfileHMM& y) = 0;
    virtual PairwiseAlignment make_alignment(const ProfileHMM& x, const ProfileHMM& y, const DPTrace& tr) = 0;
};

#endif

