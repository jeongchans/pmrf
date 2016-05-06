#ifndef _MRF_H_
#define _MRF_H_

#include <vector>
#include <string>
#include <map>
#include <unordered_map>

#include "util/common.h"
#include "seq/alphabet.h"

#define MRF_FORMAT_NUMBER 1

using std::vector;
using std::string;
using std::map;
using std::unordered_map;

typedef size_t NodeIndex;

struct EdgeIndex {
    size_t idx1;
    size_t idx2;

    EdgeIndex() {}
    EdgeIndex(const size_t& idx1, const size_t& idx2) : idx1(idx1), idx2(idx2) {}
    bool operator<(const EdgeIndex& other) const { return (idx1 < other.idx1) || (idx1 == other.idx1 && idx2 < other.idx2); }
    bool operator==(const EdgeIndex& other) const { return (idx1 == other.idx1 && idx2 == other.idx2); }
};

typedef vector<NodeIndex> NodeIndexVector;
typedef vector<EdgeIndex> EdgeIndexVector;

class MRF {

  public:
    class Element;
    class NodeElement;
    class EdgeElement;

    class Visitor {
      public:
        virtual void visit_node(NodeElement* element, const size_t& idx) = 0;
        virtual void visit_edge(EdgeElement* element, const size_t& idx1, const size_t& idx2) = 0;
    };

    class Element {
      public:
        virtual void accept(Visitor* visitor, const size_t& idx1, const size_t& idx2) = 0;
    };

    class NodeElement : public Element {
      public:
        NodeElement() : num_weight(0) {};
        NodeElement(const size_t& num_weight) : num_weight(num_weight) {};

        NodeElement& operator=(const NodeElement& rhs);

        virtual void accept(Visitor* visitor, const size_t& idx1, const size_t& idx2=(size_t) NULL) { visitor->visit_node(this, idx1); }

        // weight getter and setter
        const VectorXf& get_weight() const { return weight; }
        void set_weight(const VectorXf& w) { weight = w; }
        void set_weight(const double& w) { set_weight(VectorXf::Constant(num_weight, w)); }

      protected:
        size_t num_weight;
        VectorXf weight;
    };

    class EdgeElement : public Element {
      public:
        EdgeElement() : num_weight1(0), num_weight2(0) {};
        EdgeElement(const size_t& num_weight1, const size_t& num_weight2) : num_weight1(num_weight1), num_weight2(num_weight2) {};

        EdgeElement& operator=(const EdgeElement& rhs);

        virtual void accept(Visitor* visitor, const size_t& idx1, const size_t& idx2) { visitor->visit_edge(this, idx1, idx2); }

        // weight getter and setter
        const MatrixXf& get_weight() const { return weight; }
        void set_weight(const MatrixXf& w) { weight = w; }
        void set_weight(const double& w) { set_weight(MatrixXf::Constant(num_weight1, num_weight2, w)); }

      protected:
        size_t num_weight1, num_weight2;
        MatrixXf weight;
    };

  public:
    MRF(const size_t& length, const Alphabet& abc, const EdgeIndexVector* eidxs=NULL);
    MRF(const string& seq, const Alphabet& abc, const EdgeIndexVector* eidxs=NULL);

    void traverse(Visitor& visitor);

    // length, sequence
    string get_seq() const { return seq; }
    size_t get_length() const { return length; }
    size_t get_num_edge() const { return edges.size(); }
    size_t get_num_var() const { return abc.get_canonical_size(); }
    string get_var_symbol() const { return abc.get_canonical(); }
    const Alphabet& get_alphabet() const { return abc; }

    // Profile setter and getter
    void set_psfm(const MatrixXf& m) { psfm = m; }
    MatrixXf get_psfm() const { return psfm; }
    const VectorXf get_psfm(const size_t& idx) const { return psfm.row(idx); }

    // element getter
    NodeElement& get_node(const size_t& idx) { return nodes[idx]; }
    EdgeElement& get_edge(const size_t& idx1, const size_t& idx2) { return get_edge(EdgeIndex(idx1, idx2)); }
    EdgeElement& get_edge(const EdgeIndex& eidx) { return edges.at(eidx); }

    const NodeElement& get_node(const size_t& idx) const { return nodes[idx]; }
    const EdgeElement& get_edge(const size_t& idx1, const size_t& idx2) const 
        { return get_edge(EdgeIndex(idx1, idx2)); }
    const EdgeElement& get_edge(const EdgeIndex& eidx) const;
    EdgeIndexVector get_edge_idxs() const;

    bool has_edge(const size_t& idx1, const size_t& idx2) const { return has_edge(EdgeIndex(idx1, idx2)); }
    bool has_edge(const EdgeIndex& eidx) const;

    // modeling log
    void set_neff(const float& n) { neff = n; }
    float get_neff() const { return neff; }
    void set_fmtnum(const size_t& n) { fmtnum = n; }
    size_t get_fmtnum() const { return fmtnum; }

  private:
    const Alphabet& abc;
    string seq;
    size_t length;
    vector<NodeElement> nodes;
    map<EdgeIndex, EdgeElement> edges;
    MatrixXf psfm;  // Position-specific frequency matrix
    float neff;
    size_t fmtnum;

    void init(const EdgeIndexVector* eidxs);
};

inline bool MRF::has_edge(const EdgeIndex& eidx) const {
    if (eidx.idx1 < eidx.idx2) return edges.find(eidx) != edges.end();
    else return edges.find(EdgeIndex(eidx.idx2, eidx.idx1)) != edges.end();
}

inline const MRF::EdgeElement& MRF::get_edge(const EdgeIndex& eidx) const {
    if (eidx.idx1 < eidx.idx2) return edges.find(eidx)->second;
    else return edges.find(EdgeIndex(eidx.idx2, eidx.idx1))->second;
}

#endif
