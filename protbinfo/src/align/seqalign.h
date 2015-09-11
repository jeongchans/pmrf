#ifndef _SEQALIGN_H_
#define _SEQALIGN_H_

#include <string>

#include "dpalg.h"
#include "seq/subsmat.h"

using std::string;

class NWMatrixInitializer : public DPMatrixInitializer {
  public:
    NWMatrixInitializer(const double& gap_open, 
                        const double& gap_ext, 
                        const DPMatrix& matrix) 
    : gap_open(gap_open), 
      gap_ext(gap_ext), 
      matrix(matrix) {};

    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col);

  private:
    const double& gap_open;
    const double& gap_ext;
    const DPMatrix& matrix;
};

class NWMatrixOptimizer : public DPMatrixOptimizer {
  public:
    NWMatrixOptimizer(const SubstitutionMatrix& subsmat, 
                      const string& seq1,
                      const string& seq2,
                      const double& gap_open, 
                      const double& gap_ext, 
                      const double& baseline,
                      const DPMatrix& matrix) 
    : subsmat(subsmat), 
      seq1(seq1),
      seq2(seq2),
      gap_open(gap_open), 
      gap_ext(gap_ext), 
      baseline(baseline),
      matrix(matrix) {};

    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col);

  protected:
    virtual pair<double, DPElementPosition> find_optimum(const DPMatrix& matrix);

  private:
    const SubstitutionMatrix& subsmat;
    const string& seq1;
    const string& seq2;
    const double& gap_open;
    const double& gap_ext;
    const double& baseline;
    const DPMatrix& matrix;
};

class NWAlign : public DPAlgorithm {
  public:
    NWAlign(const SubstitutionMatrix& subsmat, const double& gap_open, const double& gap_ext, const double& baseline)
    : subsmat(subsmat), 
      gap_open(gap_open), 
      gap_ext(gap_ext), 
      baseline(baseline) {};

  private:
    const SubstitutionMatrix& subsmat;
    const double& gap_open;
    const double& gap_ext;
    const double& baseline;

    virtual void init_dp_align(const string& x, const string& y);
    virtual PairwiseAlignment make_alignment(const string& x, const string& y, const DPTrace& tr);
};

class SWMatrixInitializer : public DPMatrixInitializer {
  public:
    SWMatrixInitializer(const double& gap_open, 
                        const double& gap_ext, 
                        const DPMatrix& matrix) 
    : gap_open(gap_open), 
      gap_ext(gap_ext), 
      matrix(matrix) {};

    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col);

  private:
    const double& gap_open;
    const double& gap_ext;
    const DPMatrix& matrix;
};

class SWMatrixOptimizer : public DPMatrixOptimizer {
  public:
    SWMatrixOptimizer(const SubstitutionMatrix& subsmat, 
                      const string& seq1,
                      const string& seq2,
                      const double& gap_open, 
                      const double& gap_ext, 
                      const double& baseline,
                      const DPMatrix& matrix) 
    : subsmat(subsmat), 
      seq1(seq1),
      seq2(seq2),
      gap_open(gap_open), 
      gap_ext(gap_ext), 
      baseline(baseline),
      matrix(matrix) {};

    virtual void visit_match(MatchMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_insdel(InsDelMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_delins(DelInsMatrixElement* element, const size_t& row, const size_t& col);

  protected:
    virtual pair<double, DPElementPosition> find_optimum(const DPMatrix& matrix);

  private:
    const SubstitutionMatrix& subsmat;
    const string& seq1;
    const string& seq2;
    const double& gap_open;
    const double& gap_ext;
    const double& baseline;
    const DPMatrix& matrix;
};

class SWAlign : public DPAlgorithm {
  public:
    SWAlign(const SubstitutionMatrix& subsmat, const double& gap_open, const double& gap_ext, const double& baseline)
    : subsmat(subsmat), 
      gap_open(gap_open), 
      gap_ext(gap_ext), 
      baseline(baseline) {};

  private:
    const SubstitutionMatrix& subsmat;
    const double& gap_open;
    const double& gap_ext;
    const double& baseline;

    virtual void init_dp_align(const string& x, const string& y);
    virtual PairwiseAlignment make_alignment(const string& x, const string& y, const DPTrace& tr);
};

#endif
