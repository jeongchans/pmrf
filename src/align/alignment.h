#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <vector>
#include <map>
#include <string>

const int GAP = -1;

class Alignment {
  public:
    std::string name;
    std::string desc;
    std::string seq;            // aligned sequence 
    std::vector<int> pos_idx;   // aligned position indices

    void push_back(const char& c, const int& idx);
    int get_aligned_seq_length() const;
};

struct PairwiseAlignment {
    Alignment alignment[2];
    std::map<std::string, double> property;     // property such as score and E-value
};

#endif
