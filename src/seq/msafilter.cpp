#include "msafilter.h"

#include <algorithm>

vector<string> TerminalGapRemover::build_terminus_state(const vector<string>& msa) const {
    vector<string> st;
    for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos) {
        size_t last_idx = pos->size() - 1;
        size_t found1 = pos->find_first_not_of(abc.get_gap());
        size_t found2 = pos->find_last_not_of(abc.get_gap());
        st.push_back(string(found1, TERMI) +
                     string(found2 - found1 + 1, NON_TERMI) + 
                     string(last_idx - found2, TERMI));
    }
    return st;
}

vector<string> TerminalGapRemover::filter(const vector<string>& msa) const {
    vector<string> st = build_terminus_state(msa);;
    double gap_perc;
    int rows = msa.size();
    int cols = msa[0].size();
    int start = 0;
    int end = cols - 1;
    for (int i = 0; i < cols; ++i) {
        string s = "";
        for (vector<string>::iterator pos = st.begin(); pos != st.end(); ++pos)
            s += (*pos)[i];
        gap_perc = (double) std::count(s.begin(), s.end(), TERMI) / (double) rows;
        if (gap_perc <= max_gap_perc) {
            start = i;
            break;
        }
    }
    for (int i = cols - 1; i >= 0; --i) {
        string s = "";
        for (vector<string>::iterator pos = st.begin(); pos != st.end(); ++pos)
            s += (*pos)[i];
        gap_perc = (double) std::count(s.begin(), s.end(), TERMI) / (double) rows;
        if (gap_perc <= max_gap_perc) {
            end = i;
            break;
        }
    }
    vector<string> r;
    for (vector<string>::const_iterator pos = msa.begin(); pos != msa.end(); ++pos)
        r.push_back(pos->substr(start, end - start + 1));
    return r;
}
