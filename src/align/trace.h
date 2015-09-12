#ifndef _PROTBINFO_TRACE_H_
#define _PROTBINFO_TRACE_H_

#include <vector>
#include <string>

#include <blitz/array.h>

#include "util/common.h"
#include "util/numeric.h"
#include "seq/alphabet.h"

enum StateType { MATCH = 0, DELETE = 1, INSERT = 2, UNDEFINED = 99 };
static const size_t NUM_STATE_TYPE = 3;

typedef Float1dArray FreqVec;

class Trace {
  public:
    Trace();
    Trace(const char* st, const char* seq);
    Trace(const string& st, const string& seq);
    Trace(const string& st, const string& seq, const string& id, const string& desc);

    std::string get_id() const { return id; }
    std::string get_desc() const { return desc; }
    std::string get_seq() const;
    size_t get_length() const { return length; }
    string get_property(const string& key) const { return property.find(key)->second; }
    void set_property(const string& key, const string& value) { property.insert(make_pair(key, value)); }

    bool operator==(const Trace& rhs) const;
    bool is_passing(const StateType& type, const size_t& pos) const;
    bool has_terminal_gap(const size_t& pos) const;
    std::string get_emit(const StateType& type, const size_t& pos) const;
    FreqVec get_transit_count(const StateType& type, const size_t& pos) const;
    std::string get_MD_seq() const;
    int get_visit(const StateType& type, const size_t& pos) const;

  private:
    std::string id, desc;
    size_t length;  // length of reference sequence
    size_t first_not_delete, last_not_delete; // first and last non-deletion index
    std::vector<std::string> emit[NUM_STATE_TYPE];     // emitted symbols at each state
    map<string, string> property;

    void init_trace(const string& st, const string& seq);
};

class TraceVector {
  public:
    typedef vector<Trace>::iterator iterator;
    typedef vector<Trace>::const_iterator const_iterator;
    typedef vector<Trace>::reverse_iterator reverse_iterator;

    TraceVector() {};

    Trace& operator[](const size_t& n) { return data[n]; }
    const Trace& operator[](const size_t& n) const { return data[n]; }
    size_t size() const { return data.size(); }
    iterator begin() { return data.begin(); }
    const_iterator begin() const { return data.begin(); }
    iterator end() { return data.end(); }
    const_iterator end() const { return data.end(); }
    void push_back(const Trace& x) { data.push_back(x); }
    bool empty() const { return data.empty(); }
    reverse_iterator rbegin() { return data.rbegin(); }

    TraceVector subset_passing(const StateType& type, const size_t& idx, const bool& omit_termi_gap=true) const;
    vector<string> get_MD_seq_vec() const;

    string ss_dssp;
    string ss_pred;
    string ss_conf;

  private:
    vector<Trace> data;
};

enum MSAFormat { AFASTA, A3M };

class TraceImporter {
  public:
    TraceImporter(const Alphabet& abc) : abc(abc) {};

    std::pair<size_t, TraceVector> import(std::istream& is, const MSAFormat& fmt);

  private:
    const Alphabet& abc;

    void afa_to_trace(std::istream& is, TraceVector& traces);
    void a3m_to_trace(std::istream& is, TraceVector& traces);
};

#endif
