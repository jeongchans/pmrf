#include "bgfreq.h"

VectorXf BgFreq::get_array(const Alphabet& abc) const {
    size_t n = abc.get_canonical_size();
    VectorXf v(n);
    int i;
    for (std::map<char, double>::const_iterator pos = freq.begin(); pos != freq.end(); ++pos) {
        i = abc.get_idx(pos->first);
        v(i) = pos->second;
    }
    return v;
}

RobinsonBgFreq::RobinsonBgFreq() {
    freq['A'] = 0.07805;
    freq['C'] = 0.01925;
    freq['D'] = 0.05364;
    freq['E'] = 0.06295;
    freq['F'] = 0.03856;
    freq['G'] = 0.07377;
    freq['H'] = 0.02199;
    freq['I'] = 0.05142;
    freq['K'] = 0.05744;
    freq['L'] = 0.09019;
    freq['M'] = 0.02243;
    freq['N'] = 0.04487;
    freq['P'] = 0.05203;
    freq['Q'] = 0.04264;
    freq['R'] = 0.05129;
    freq['S'] = 0.07120;
    freq['T'] = 0.05841;
    freq['V'] = 0.06441;
    freq['W'] = 0.01330;
    freq['Y'] = 0.03216;
}

AltschulBgFreq::AltschulBgFreq() {
    freq['A'] = 0.08100;
    freq['C'] = 0.01500;
    freq['D'] = 0.05400;
    freq['E'] = 0.06100;
    freq['F'] = 0.04000;
    freq['G'] = 0.06800;
    freq['H'] = 0.02200;
    freq['I'] = 0.05700;
    freq['K'] = 0.05600;
    freq['L'] = 0.09300;
    freq['M'] = 0.02500;
    freq['N'] = 0.04500;
    freq['P'] = 0.04900;
    freq['Q'] = 0.03900;
    freq['R'] = 0.05700;
    freq['S'] = 0.06800;
    freq['T'] = 0.05800;
    freq['V'] = 0.06700;
    freq['W'] = 0.01300;
    freq['Y'] = 0.03200;
}
