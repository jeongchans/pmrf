#ifndef _SEQ_SUBSMAT_H_
#define _SEQ_SUBSMAT_H_

#include "util/common.h"
#include "util/numeric.h"
#include "alphabet.h"
#include "rootfinder.h"

class SubstitutionMatrix {
  public:
    SubstitutionMatrix() {};
    SubstitutionMatrix(const std::map<SymbolPair, double>& score) : score(score) {};

    MatrixXf get_array(const Alphabet& abc) const;
    double get_value(const char& x, const char& y) const {
        std::map<SymbolPair, double>::const_iterator pos;
        pos = score.find(SymbolPair(x, y));
        if (pos != score.end()) return pos->second;
        else return score.find(SymbolPair(y, x))->second;
    }

  protected:
    std::map<SymbolPair, double> score;
};

class BLOSUM62Matrix : public SubstitutionMatrix {
  public:
    BLOSUM62Matrix();
};

class GonnetMatrix : public SubstitutionMatrix {
  public:
    GonnetMatrix();
};

class TargetProbEstimatorGivenBG {

  private:
    class TargetProbFunction : public TargetFunction {
      public:
        TargetProbFunction(const VectorXf& bgfreq, const MatrixXf& scoremat) {
            this->mfreq = (bgfreq * bgfreq.transpose()).cast<double>();
            this->scoremat = scoremat.cast<double>();
        }

        virtual double fx(const double& x) const {
            return mfreq.cwiseProduct((x * scoremat).unaryExpr(&exp)).sum() - 1.;
        }

        virtual pair<double, double> fx_dfx(const double& x) const {
            return make_pair(fx(x), mfreq.cwiseProduct((x * scoremat).unaryExpr(&exp)).cwiseProduct(scoremat).sum());
        }

      private:
        MatrixXd mfreq;        // marginal frequency
        MatrixXd scoremat;     // substitution score
    };

  public:
    TargetProbEstimatorGivenBG(const VectorXf& bgfreq) {
        this->bgfreq = bgfreq;
    }

    pair<double, MatrixXf> probify(const MatrixXf& scoremat);

  private:
    VectorXf bgfreq;
};

static const BLOSUM62Matrix BLOSUM62_MATRIX;

#endif
