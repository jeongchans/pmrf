#include "mrf.h"

/* NodeElement */
MRF::NodeElement& MRF::NodeElement::operator=(const NodeElement& rhs) {
    num_weight = rhs.num_weight;
    set_weight(rhs.weight);
    return *this;
}

/* EdgeElement */
MRF::EdgeElement& MRF::EdgeElement::operator=(const EdgeElement& rhs) {
    num_weight1 = rhs.num_weight1;
    num_weight2 = rhs.num_weight2;
    set_weight(rhs.weight);
    return *this;
}

/* MRF */
MRF::MRF(const size_t& length, const Alphabet& abc, const EdgeIndexVector* eidxs) 
: abc(abc), length(length), fmtnum(MRF_FORMAT_NUMBER) {
    seq = "";
    char c = abc.get_unknown()[0];
    for (size_t i = 0; i < length; ++i) seq += c;
    init(eidxs);
}

MRF::MRF(const string& seq, const Alphabet& abc, const EdgeIndexVector* eidxs) 
: abc(abc), seq(seq), fmtnum(MRF_FORMAT_NUMBER) {
    length = seq.size();
    init(eidxs);
}

void MRF::traverse(MRF::Visitor& visitor) {
    for (size_t idx = 0; idx < length; ++idx)
        nodes[idx].accept(&visitor, idx);
    for (map<EdgeIndex, EdgeElement>::iterator pos = edges.begin(); pos != edges.end(); ++pos)
        (pos->second).accept(&visitor, pos->first.idx1, pos->first.idx2);
}

EdgeIndexVector MRF::get_edge_idxs() const {
    EdgeIndexVector idxs;
    for (map<EdgeIndex, EdgeElement>::const_iterator pos = edges.begin(); pos != edges.end(); ++pos)
        idxs.push_back(EdgeIndex(pos->first.idx1, pos->first.idx2));
    return idxs;
}

void MRF::init(const EdgeIndexVector* eidxs) {
    nodes.clear();
    edges.clear();
    size_t n = abc.get_canonical_size();
    for (size_t i = 0; i < length; ++i) nodes.push_back(NodeElement(n));
    if (eidxs != NULL) {
        for (EdgeIndexVector::const_iterator pos = eidxs->begin(); pos != eidxs->end(); ++pos)
            edges[*pos] = EdgeElement(n, n);
    } else {
        for (size_t i = 0; i < length; ++i)
            for (size_t j = i + 1; j < length; ++j)
                edges[EdgeIndex(i, j)] = EdgeElement(n, n);
    }
    psfm = MatrixXf::Zero(length, n);
}
