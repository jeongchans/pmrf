#include "subsmat.h"

MatrixXf SubstitutionMatrix::get_array(const Alphabet& abc) const {
    size_t n = abc.get_canonical_size();
    MatrixXf m(n, n);
    int i, j;
    for (std::map<SymbolPair, double>::const_iterator pos = score.begin(); pos != score.end(); ++pos) {
        i = abc.get_idx((pos->first).first);
        j = abc.get_idx((pos->first).second);
        m(i, j) = pos->second;
        m(j, i) = pos->second;
    }
    return m;
}

BLOSUM62Matrix::BLOSUM62Matrix() {
    std::string s = "ACDEFGHIKLMNPQRSTVWY";
    MatrixXf m(20, 20);
    //     A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
    m <<   4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2,  // A
           0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2,  // C
          -2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3,  // D
          -1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2,  // E
          -2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3,  // F
           0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3,  // G
          -2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2,  // H
          -1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1,  // I
          -1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2,  // K
          -1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1,  // L
          -1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1,  // M
          -2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2,  // N
          -1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3,  // P
          -1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1,  // Q
          -1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2,  // R
           1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2,  // S
           0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2,  // T
           0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1,  // V
          -3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2,  // W
          -2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7;  // Y
    for (int i = 0; i < 20; ++i)
        for (int j = i; j < 20; ++j)
            score.insert(make_pair(SymbolPair(s[i], s[j]), m(i, j)));
}

GonnetMatrix::GonnetMatrix() {
    string s = "ARNDCQEGHILKMFPSTWYVBZX*";
    int n = s.size();
    MatrixXf m(n, n);
    //    A     R     N     D     C     Q     E     G     H     I     L     K     M     F     P     S     T     W     Y     V     B     Z     X     *     
    m <<  2.4, -0.6, -0.3, -0.3,  0.5, -0.2,  0.0,  0.5, -0.8, -0.8, -1.2, -0.4, -0.7, -2.3,  0.3,  1.1,  0.6, -3.6, -2.2,  0.1, -0.5,  0.0, -0.4, -8.0,  // A
         -0.6,  4.7,  0.3, -0.3, -2.2,  1.5,  0.4, -1.0,  0.6, -2.4, -2.2,  2.7, -1.7, -3.2, -0.9, -0.2, -0.2, -1.6, -1.8, -2.0, -0.2,  1.0, -0.5, -8.0,  // R
         -0.3,  0.3,  3.8,  2.2, -1.8,  0.7,  0.9,  0.4,  1.2, -2.8, -3.0,  0.8, -2.2, -3.1, -0.9,  0.9,  0.5, -3.6, -1.4, -2.2,  2.8,  0.9, -0.2, -8.0,  // N
         -0.3, -0.3,  2.2,  4.7, -3.2,  0.9,  2.7,  0.1,  0.4, -3.8, -4.0,  0.5, -3.0, -4.5, -0.7,  0.5,  0.0, -5.2, -2.8, -2.9,  3.5,  2.2, -0.5, -8.0,  // D
          0.5, -2.2, -1.8, -3.2, 11.5, -2.4, -3.0, -2.0, -1.3, -1.1, -1.5, -2.8, -0.9, -0.8, -3.1,  0.1, -0.5, -1.0, -0.5,  0.0, -2.7, -2.6, -1.1, -8.0,  // C
         -0.2,  1.5,  0.7,  0.9, -2.4,  2.7,  1.7, -1.0,  1.2, -1.9, -1.6,  1.5, -1.0, -2.6, -0.2,  0.2,  0.0, -2.7, -1.7, -1.5,  0.6,  2.3, -0.1, -8.0,  // Q
          0.0,  0.4,  0.9,  2.7, -3.0,  1.7,  3.6, -0.8,  0.4, -2.7, -2.8,  1.2, -2.0, -3.9, -0.5,  0.2, -0.1, -4.3, -2.7, -1.9,  1.8,  3.1, -0.4, -8.0,  // E
          0.5, -1.0,  0.4,  0.1, -2.0, -1.0, -0.8,  6.6, -1.4, -4.5, -4.4, -1.1, -3.5, -5.2, -1.6,  0.4, -1.1, -4.0, -4.0, -3.3,  0.0, -0.8, -1.5, -8.0,  // G
         -0.8,  0.6,  1.2,  0.4, -1.3,  1.2,  0.4, -1.4,  6.0, -2.2, -1.9,  0.6, -1.3, -0.1, -1.1, -0.2, -0.3, -0.8,  2.2, -2.0,  0.6,  0.9, -0.2, -8.0,  // H
         -0.8, -2.4, -2.8, -3.8, -1.1, -1.9, -2.7, -4.5, -2.2,  4.0,  2.8, -2.1,  2.5,  1.0, -2.6, -1.8, -0.6, -1.8, -0.7,  3.1, -3.5, -2.2, -0.5, -8.0,  // I
         -1.2, -2.2, -3.0, -4.0, -1.5, -1.6, -2.8, -4.4, -1.9,  2.8,  4.0, -2.1,  2.8,  2.0, -2.3, -2.1, -1.3, -0.7,  0.0,  1.8, -3.7, -2.2, -1.0, -8.0,  // L
         -0.4,  2.7,  0.8,  0.5, -2.8,  1.5,  1.2, -1.1,  0.6, -2.1, -2.1,  3.2, -1.4, -3.3, -0.6,  0.1,  0.1, -3.5, -2.1, -1.7,  0.4,  1.5, -0.3, -8.0,  // K
         -0.7, -1.7, -2.2, -3.0, -0.9, -1.0, -2.0, -3.5, -1.3,  2.5,  2.8, -1.4,  4.3,  1.6, -2.4, -1.4, -0.6, -1.0, -0.2,  1.6, -2.8, -1.4, -0.1, -8.0,  // M
         -2.3, -3.2, -3.1, -4.5, -0.8, -2.6, -3.9, -5.2, -0.1,  1.0,  2.0, -3.3,  1.6,  7.0, -3.8, -2.8, -2.2,  3.6,  5.1,  0.1, -4.0, -3.2, -0.8, -8.0,  // F
          0.3, -0.9, -0.9, -0.7, -3.1, -0.2, -0.5, -1.6, -1.1, -2.6, -2.3, -0.6, -2.4, -3.8,  7.6,  0.4,  0.1, -5.0, -3.1, -1.8, -1.0, -0.2, -1.1, -8.0,  // P
          1.1, -0.2,  0.9,  0.5,  0.1,  0.2,  0.2,  0.4, -0.2, -1.8, -2.1,  0.1, -1.4, -2.8,  0.4,  2.2,  1.5, -3.3, -1.9, -1.0,  0.5,  0.3, -0.2, -8.0,  // S
          0.6, -0.2,  0.5,  0.0, -0.5,  0.0, -0.1, -1.1, -0.3, -0.6, -1.3,  0.1, -0.6, -2.2,  0.1,  1.5,  2.5, -3.5, -1.9,  0.0,  0.0,  0.1, -0.2, -8.0,  // T
         -3.6, -1.6, -3.6, -5.2, -1.0, -2.7, -4.3, -4.0, -0.8, -1.8, -0.7, -3.5, -1.0,  3.6, -5.0, -3.3, -3.5, 14.2,  4.1, -2.6, -4.6, -3.5, -1.6, -8.0,  // W
         -2.2, -1.8, -1.4, -2.8, -0.5, -1.7, -2.7, -4.0,  2.2, -0.7,  0.0, -2.1, -0.2,  5.1, -3.1, -1.9, -1.9,  4.1,  7.8, -1.1, -2.3, -2.1, -0.7, -8.0,  // Y
          0.1, -2.0, -2.2, -2.9,  0.0, -1.5, -1.9, -3.3, -2.0,  3.1,  1.8, -1.7,  1.6,  0.1, -1.8, -1.0,  0.0, -2.6, -1.1,  3.4, -2.8, -1.6, -0.5, -8.0,  // V
         -0.5, -0.2,  2.8,  3.5, -2.7,  0.6,  1.8,  0.0,  0.6, -3.5, -3.7,  0.4, -2.8, -4.0, -1.0,  0.5,  0.0, -4.6, -2.3, -2.8,  3.0,  1.5, -2.6, -8.0,  // B
          0.0,  1.0,  0.9,  2.2, -2.6,  2.3,  3.1, -0.8,  0.9, -2.2, -2.2,  1.5, -1.4, -3.2, -0.2,  0.3,  0.1, -3.5, -2.1, -1.6,  1.5,  2.9, -2.2, -8.0,  // Z
         -0.4, -0.5, -0.2, -0.5, -1.1, -0.1, -0.4, -1.5, -0.2, -0.5, -1.0, -0.3, -0.1, -0.8, -1.1, -0.2, -0.2, -1.6, -0.7, -0.5, -2.6, -2.2, -0.3, -8.0,  // X
         -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0, -8.0,  1.0;  // *
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j)
            score.insert(make_pair(SymbolPair(s[i], s[j]), m(i, j)));
}

pair<double, MatrixXf> TargetProbEstimatorGivenBG::probify(const MatrixXf& scoremat) {
    // initial guest at lambda
    TargetProbFunction target_func(bgfreq, scoremat);
    double lamb;
    double fx;
    for (lamb = 1. / scoremat.maxCoeff(); lamb < 50.; lamb *= 2.) {
        fx = target_func.fx(lamb);
        if (fx > 0) break;
    }
    // use Newton/Raphson method to find the optimal lambda
    NewtonRaphsonRootFinder solver;
    lamb = solver.find_root(target_func, lamb);
    MatrixXf prob(scoremat.rows(), scoremat.cols());
    prob = (bgfreq * bgfreq.transpose()).cwiseProduct((lamb * scoremat).unaryExpr(&exp));
    return make_pair(lamb, prob);
}
