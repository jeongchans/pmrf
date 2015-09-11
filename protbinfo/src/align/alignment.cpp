#include "alignment.h"

void Alignment::push_back(const char& c, const int& idx) {
    seq += c;
    pos_idx.push_back(idx);
}

int Alignment::get_aligned_seq_length() const { 
    return *pos_idx.rbegin() - *pos_idx.begin() + 1;
}
