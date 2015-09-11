#ifndef _DPBASE_H_
#define _DPBASE_H_

#include <deque>

using std::deque;
using std::pair;

class DPMatrixElementBase;

struct DPElementPosition {
    const DPMatrixElementBase* element;
    int row;    // row index
    int col;    // column index
};

class DPMatrixElementBase {
  public:
    double score;
    DPElementPosition from;
};

struct DPTrace {
    deque<pair<int, int> > res_trace;
    deque<DPElementPosition> element_trace;
};

#endif
