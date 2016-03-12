#ifndef _MRF_H_
#define _MRF_H_

#include <vector>
#include <string>
#include <map>

#include <gtest/gtest_prod.h>

#include "util/common.h"
#include "seq/alphabet.h"

using std::vector;
using std::string;
using std::map;

struct EdgeIndex {
    size_t idx1;
    size_t idx2;

    EdgeIndex(const size_t& idx1, const size_t& idx2) : idx1(idx1), idx2(idx2) {}
    bool operator<(const EdgeIndex& other) const { return (idx1 < other.idx1) || (idx1 == other.idx1 && idx2 < other.idx2); }
    bool operator==(const EdgeIndex& other){ return (idx1 == other.idx1 && idx2 == other.idx2) || (idx1 == other.idx2 && idx2 == other.idx1); }
};

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

        virtual void accept(Visitor* visitor, const size_t& idx1, const size_t& idx2);

        // weight getter and setter
        const Float1dArray& get_weight() const;
        void set_weight(const Float1dArray& w);
        void set_weight(const double& w);

      protected:
        size_t num_weight;
        Float1dArray weight;
    };

    class EdgeElement : public Element {
      public:
        EdgeElement() : num_weight1(0), num_weight2(0) {};
        EdgeElement(const size_t& num_weight1, const size_t& num_weight2) : num_weight1(num_weight1), num_weight2(num_weight2) {};

        EdgeElement& operator=(const EdgeElement& rhs);

        virtual void accept(Visitor* visitor, const size_t& idx1, const size_t& idx2=(size_t) NULL);

        // weight getter and setter
        const Float2dArray& get_weight() const;
        void set_weight(const Float2dArray& w);
        void set_weight(const double& w);

      protected:
        size_t num_weight1, num_weight2;
        Float2dArray weight;
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
    void set_psfm(const Float2dArray& m) { psfm = m; }
    Float2dArray get_psfm() const { return psfm; }
    const Float1dArray get_psfm(const size_t& idx) const { return psfm(idx, ALL); }

    // element getter
    NodeElement& get_node(const size_t& idx) { return nodes[idx]; }
    EdgeElement& get_edge(const size_t& idx1, const size_t& idx2) { return edges[EdgeIndex(idx1, idx2)]; }

    const NodeElement& get_node(const size_t& idx) const
        { return nodes[idx]; }
    const EdgeElement& get_edge(const size_t& idx1, const size_t& idx2) const
        { return edges.find(EdgeIndex(idx1, idx2))->second; }
    EdgeIndexVector get_edge_idxs() const;

    bool has_edge(const size_t& idx1, const size_t& idx2) const
        { return edges.find(EdgeIndex(idx1, idx2)) != edges.end(); }

  private:
    const Alphabet& abc;
    string seq;
    size_t length;
    vector<NodeElement> nodes;
    map<EdgeIndex, EdgeElement> edges;
    Float2dArray psfm;  // Position-specific frequency matrix

    void init(const EdgeIndexVector* eidxs);
};

#endif
