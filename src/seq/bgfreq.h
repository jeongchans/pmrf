#ifndef _SEQ_BGFREQ_H_
#define _SEQ_BGFREQ_H_

#include "util/common.h"
#include "alphabet.h"

class BgFreq {
  public:
    BgFreq() {};
    BgFreq(const std::map<char, double>& freq) : freq(freq) {};

    VectorXf get_array(const Alphabet& abc) const;
    double get_value(const char& x) const { return freq.find(x)->second; }

  protected:
    std::map<char, double> freq;
};

// Robinson and Robinson (default in PSI-BLAST)
// PNAS-88-8880(1991)
class RobinsonBgFreq : public BgFreq {
  public:
    RobinsonBgFreq();
};

// Altschul's background frequency
class AltschulBgFreq : public BgFreq {
  public:
    AltschulBgFreq();
};

static const RobinsonBgFreq ROBINSON_BGFREQ;

#endif
