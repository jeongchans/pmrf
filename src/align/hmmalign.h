#ifndef _HMMALIGN_H_
#define _HMMALIGN_H_

#include "hhdpalg.h"
#include "seq/bgfreq.h"

/**
 * Local HMM-HMM alignment
 */

class LocalHHMatrixInitializer : public HHDPMatrixInitializer {
  public:
    LocalHHMatrixInitializer(const ProfileHMM& hmm1,
                             const ProfileHMM& hmm2,
                             const HHDPMatrix& matrix) 
    : hmm1(hmm1),
      hmm2(hmm2),
      matrix(matrix) {};

    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col);

  private:
    const ProfileHMM& hmm1;
    const ProfileHMM& hmm2;
    const HHDPMatrix& matrix;
};

class LocalHHMatrixOptimizer : public HHDPMatrixOptimizer {
  public:
    LocalHHMatrixOptimizer(const Float1dArray& bgfreq,
                           const ProfileHMM& hmm1,
                           const ProfileHMM& hmm2,
                           const double& baseline,
                           const HHDPMatrix& matrix) 
    : bgfreq(bgfreq),
      hmm1(hmm1),
      hmm2(hmm2),
      baseline(baseline),
      matrix(matrix) {};

    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col);

  protected:
    virtual pair<double, DPElementPosition> find_optimum(const HHDPMatrix& matrix);

    double calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to);
    double calc_column_score(const HMMMatchState& match1, const HMMMatchState& match2);

  private:
    const Float1dArray bgfreq;
    const ProfileHMM& hmm1;
    const ProfileHMM& hmm2;
    const double& baseline;
    const HHDPMatrix& matrix;
};

class LocalHHAlign : public HHDPAlgorithm {
  public:
    LocalHHAlign(const Float1dArray& bgfreq, const double& baseline)
    : bgfreq(bgfreq), 
      baseline(baseline) {};

  private:
    const Float1dArray bgfreq;
    const double& baseline;

    virtual void init_dp_align(const ProfileHMM& x, const ProfileHMM& y);
    virtual PairwiseAlignment make_alignment(const ProfileHMM& hmm1, const ProfileHMM& hmm2, const DPTrace& tr);
};

/**
 * Global HMM-HMM alignment
 */

class GlobalHHMatrixInitializer : public HHDPMatrixInitializer {
  public:
    GlobalHHMatrixInitializer(const ProfileHMM& hmm1,
                              const ProfileHMM& hmm2,
                              const HHDPMatrix& matrix) 
    : hmm1(hmm1),
      hmm2(hmm2),
      matrix(matrix) {};

    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col);

  protected:
    double calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to);

  private:
    const ProfileHMM& hmm1;
    const ProfileHMM& hmm2;
    const HHDPMatrix& matrix;
};

class GlobalHHMatrixOptimizer : public HHDPMatrixOptimizer {
  public:
    GlobalHHMatrixOptimizer(const Float1dArray& bgfreq,
                            const ProfileHMM& hmm1,
                            const ProfileHMM& hmm2,
                            const double& baseline,
                            const HHDPMatrix& matrix) 
    : bgfreq(bgfreq),
      hmm1(hmm1),
      hmm2(hmm2),
      baseline(baseline),
      matrix(matrix) {};

    virtual void visit_mm(MMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_mi(MIMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_im(IMMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_dg(DGMatrixElement* element, const size_t& row, const size_t& col);
    virtual void visit_gd(GDMatrixElement* element, const size_t& row, const size_t& col);

  protected:
    virtual pair<double, DPElementPosition> find_optimum(const HHDPMatrix& matrix);

    double calc_transit_score(const ProfileHMM& hmm, const int& idx, const StateType& from, const StateType& to);
    double calc_column_score(const HMMMatchState& match1, const HMMMatchState& match2);

    const Float1dArray bgfreq;
    const ProfileHMM& hmm1;
    const ProfileHMM& hmm2;
    const double& baseline;
    const HHDPMatrix& matrix;
};

class GlobalHHAlign : public HHDPAlgorithm {
  public:
    GlobalHHAlign(const Float1dArray& bgfreq, const double& baseline)
    : bgfreq(bgfreq), 
      baseline(baseline) {};

  protected:
    const Float1dArray bgfreq;
    const double& baseline;

  private:
    virtual void init_dp_align(const ProfileHMM& x, const ProfileHMM& y);
    virtual PairwiseAlignment make_alignment(const ProfileHMM& hmm1, const ProfileHMM& hmm2, const DPTrace& tr);
};

#endif
