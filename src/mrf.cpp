#include "mrf.h"

inline void resize_and_fill(Float1dArray& v, const size_t& size, const double& value) {
    v.resize(size);
    v = value;
}

inline void resize_and_fill(Float2dArray& m, const size_t& row_size, const size_t& col_size, const double& value) {
    m.resize(row_size, col_size);
    m = value;
}

/* NodeElement */
MRF::NodeElement& MRF::NodeElement::operator=(const NodeElement& rhs) {
    num_weight = rhs.num_weight;
    set_weight(rhs.weight);
    return *this;
}

void MRF::NodeElement::accept(Visitor* visitor, const size_t& idx1, const size_t& idx2=(size_t)NULL) {
    visitor->visit_node(this, idx1);
}

const Float1dArray& MRF::NodeElement::get_weight() const {
    return weight;
}

void MRF::NodeElement::set_weight(const Float1dArray& w) {
    weight.resize(w.size());
    weight = w;
}

void MRF::NodeElement::set_weight(const double& w) {
    resize_and_fill(weight, num_weight, w);
}

/* EdgeElement */
MRF::EdgeElement& MRF::EdgeElement::operator=(const EdgeElement& rhs) {
    num_weight1 = rhs.num_weight1;
    num_weight2 = rhs.num_weight2;
    set_weight(rhs.weight);
    return *this;
}

void MRF::EdgeElement::accept(Visitor* visitor, const size_t& idx1, const size_t& idx2) {
    visitor->visit_edge(this, idx1, idx2);
}

const Float2dArray& MRF::EdgeElement::get_weight() const {
    return weight;
}

void MRF::EdgeElement::set_weight(const Float2dArray& w) {
    weight.resize(w.shape());
    weight = w;
}

void MRF::EdgeElement::set_weight(const double& w) {
    resize_and_fill(weight, num_weight1, num_weight2, w);
}

/* MRF */
MRF::MRF(const size_t& length, const Alphabet& abc, const EdgeIndexVector* eidxs) : length(length), abc(abc) {
    seq = "";
    char c = abc.get_unknown()[0];
    for (size_t i = 0; i < length; ++i) seq += c;
    init(eidxs);
}

MRF::MRF(const string& seq, const Alphabet& abc, const EdgeIndexVector* eidxs) : seq(seq), abc(abc) {
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
    resize_and_fill(psfm, length, n, 0);
}
