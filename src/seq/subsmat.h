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

    Float2dArray get_array(const Alphabet& abc) const;
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
        TargetProbFunction(const Float1dArray& bgfreq, const Float2dArray& scoremat) {
            this->mfreq.resize(bgfreq.size(), bgfreq.size());
            this->mfreq = outer(bgfreq, bgfreq);
            this->scoremat.resize(scoremat.shape());
            this->scoremat = scoremat;
        }

        virtual double fx(const double& x) const {
            return sum(mfreq * exp(x * scoremat)) - 1.;
        }
        virtual pair<double, double> fx_dfx(const double& x) const {
            return make_pair(fx(x), sum(mfreq * exp(x * scoremat) * scoremat));
        }

      private:
        Float2dArray mfreq;        // marginal frequency
        Float2dArray scoremat;     // substitution score
    };

  public:
    TargetProbEstimatorGivenBG(const Float1dArray& bgfreq) {
        this->bgfreq.resize(bgfreq.shape());
        this->bgfreq = bgfreq;
    }

    pair<double, Float2dArray> probify(const Float2dArray& scoremat);

  private:
    Float1dArray bgfreq;
};

static const BLOSUM62Matrix BLOSUM62_MATRIX;

#endif
