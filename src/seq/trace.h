#ifndef _ALIGN_TRACE_H_
#define _ALIGN_TRACE_H_

#include "util/common.h"
#include "seq/alphabet.h"

class Trace {
  public:

//    enum StateType {
//        MATCH = 0,
//        INSERT = 1,
//        GAPOPEN = 2,
//        GAPEXT = 3,
//        UNALIGNED = 4,
//        UNDEFINED = 99
//    };
//    static const size_t NUM_STATE_TYPE = 5;
//    static const char STATE_SYMBOL[5];

    enum StateType {
        MATCH = 0,
        INSERT = 1,
        DELETE = 2,
        UNDEFINED = 99
    };
    static const char STATE_SYMBOL[3];
    static const size_t NUM_STATE_TYPE = 3;

    class State {
      public:
        State(const StateType& st, const char& aa) : st(st), aa(aa) {};

        bool operator==(const State& rhs) const { return st == rhs.st && aa == rhs.aa; }
        bool operator!=(const State& rhs) const { return !(*this == rhs); }

        StateType st;
        char aa;
    };

    Trace() {};
    Trace(const string& ststr, const string& seq);
    Trace(const string& ststr, const string& seq, const string& id, const string& desc);

    std::string get_seq() const;
    size_t get_length() const { return get_seq().size(); }

    bool operator==(const Trace& rhs) const;
    std::string get_matched_aseq() const;

  private:
    std::string id, desc;
    std::vector<State> st;
    std::vector<size_t> aidx;       // converting reference index to trace index

    void init(const string& ststr, const string& seq, const string& id, const string& desc);
    StateType char_to_state(const char& x);
};

class TraceVector {
  public:
    typedef vector<Trace>::const_iterator const_iterator;

    Trace& operator[](const size_t& n) { return data[n]; }
    const Trace& operator[](const size_t& n) const { return data[n]; }
    size_t size() const { return data.size(); }
    const_iterator begin() const { return data.begin(); }
    const_iterator end() const { return data.end(); }
    void push_back(const Trace& x) { data.push_back(x); }
    bool empty() const { return data.empty(); }

    vector<string> get_matched_aseq_vec() const;

  private:
    std::vector<Trace> data;
};

enum MSAFormat { AFASTA, A3M };

class TraceImporter {
  public:
    TraceImporter(const Alphabet& abc) : abc(abc) {};

    std::pair<size_t, TraceVector> import(std::istream& is, const MSAFormat& fmt);

  private:
    const Alphabet& abc;

    void parse_state(const std::string& aseq, std::string& st, std::string& seq);
    bool is_valid_state(const std::string& st) const;
    void afa_to_trace(std::istream& is, TraceVector& traces);
    void a3m_to_trace(std::istream& is, TraceVector& traces);
};

#endif
