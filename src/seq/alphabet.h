#ifndef _ALPHABET_H_
#define _ALPHABET_H_

#include <blitz/array.h>

#include <string>
#include <map>

#include "common/numeric.h"

class Alphabet {

    class SymbolGroup {
      public:
        bool has_symbol(const char& x) const;
        void set_member(const char* symbols);
        std::string get_member() const { return member; }
        size_t get_uniq_size() const { return uniq_size; }
        void add_dup_symbol(const char& x) { member += x; }

      private:
        std::string member;
        size_t uniq_size;
    };

    typedef std::map<char, size_t> SymbolIdxMap;
    typedef std::multimap<char, char> DegeneracyMap;

  public:
    Alphabet();
    Alphabet(const char* canonical, const char* gap, const char* degenerate,
             const char* unknown, const char* none, const char* missing,
             const bool& nocase = false);

    std::string get_canonical() const { return sym_grp[CANONICAL].get_member(); }
    std::string get_gap() const { return sym_grp[GAP].get_member(); }
    std::string get_unknown() const { return sym_grp[UNKNOWN].get_member(); }
    bool is_canonical(const char& x) const;
    bool is_gap(const char& x) const;
    bool is_degenerate(const char& x) const;
    bool is_unknown(const char& x) const;
    bool is_none(const char& x) const;
    bool is_missing(const char& x) const;
    bool is_valid(const char& x) const;
    int get_idx(const char& x) const;
    size_t get_canonical_size() const;
    size_t get_valid_size() const { return get_valid_symbol().size(); }
    Float1dArray get_count(const char& x) const;
    std::string get_degeneracy(const char& x, FloatType* w=NULL) const;
    void set_degeneracy(const char& degen_ch, const char& canoni_ch);

  private:
    static const size_t CANONICAL = 0;
    static const size_t GAP = 1;
    static const size_t DEGENERATE = 2;
    static const size_t UNKNOWN = 3;
    static const size_t NONE = 4;
    static const size_t MISSING = 5;
    static const size_t NUM_GROUP = 6;

    SymbolGroup sym_grp[NUM_GROUP];
    SymbolIdxMap sym_idx;
    DegeneracyMap degeneracy;
    bool nocase;

    void update_sym_idx();
    std::string get_valid_symbol() const;

};

class AminoAcid : public Alphabet {
  public:
    AminoAcid();
};

typedef std::pair<char, char> SymbolPair;

#endif
