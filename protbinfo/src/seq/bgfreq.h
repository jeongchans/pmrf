#ifndef _BGFREQ_H_
#define _BGFREQ_H_

#include "common/common.h"
#include "common/numeric.h"
#include "alphabet.h"

class BgFreq {
  public:
    BgFreq() {};
    BgFreq(const std::map<char, FloatType>& freq) : freq(freq) {};

    Float1dArray get_array(const Alphabet& abc) const;
    FloatType get_value(const char& x) const { return freq.find(x)->second; }

  protected:
    std::map<char, FloatType> freq;
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

#endif
