#include "alphabet.h"

#include <string>
#include <algorithm>
#include <vector>
#include <iostream>

#include "util/common.h"

/**
   @class Alphabet::SymbolGroup
 */

bool Alphabet::SymbolGroup::has_symbol(const char& x) const {
    std::string::size_type idx = this->member.find(x);
    if (idx != std::string::npos) return true;
    else return false;
}

void Alphabet::SymbolGroup::set_member(const char* symbols) {
    this->member = symbols;
    this->uniq_size = this->member.size();
}

/**
   @class Alphabet
 */

Alphabet::Alphabet(const char* canonical, const char* gap, const char* degenerate,
                   const char* unknown, const char* none, const char* missing,
                   const bool& nocase, const bool& gapres) {
    sym_grp[CANONICAL].set_member(canonical);
    sym_grp[GAP].set_member(gap);
    sym_grp[DEGENERATE].set_member(degenerate);
    sym_grp[UNKNOWN].set_member(unknown);
    sym_grp[NONE].set_member(none);
    sym_grp[MISSING].set_member(missing);
    update_sym_idx();
    this->nocase = nocase;
    if (this->nocase) {
        std::string s;
        std::string::iterator pos;
        char x;
        for (size_t g_idx = 0; g_idx < NUM_GROUP; ++g_idx) {
            s = sym_grp[g_idx].get_member();
            for (pos = s.begin(); pos != s.end(); ++pos) {
                if (islower(*pos)) x = toupper(*pos);
                else if (isupper(*pos)) x = tolower(*pos);
                else continue;
                sym_grp[g_idx].add_dup_symbol(x);
                sym_idx[x] = get_idx(*pos);
            }
        }
    }
    this->gapres = gapres;
}

std::string Alphabet::get_canonical(const bool& gapres) const {
    if (gapres) return sym_grp[CANONICAL].get_member() + get_gap();
    else return sym_grp[CANONICAL].get_member();
}

size_t Alphabet::get_canonical_size(const bool& gapres) const {
    if (gapres) return sym_grp[CANONICAL].get_uniq_size() + get_gap_size();
    else return sym_grp[CANONICAL].get_uniq_size();
}

void Alphabet::set_degeneracy(const char& degen_ch, const char& canoni_ch) {
    degeneracy.insert(std::make_pair(degen_ch, canoni_ch));
}

bool Alphabet::is_canonical(const char& x, const bool& gapres) const {
    if (gapres) return sym_grp[CANONICAL].has_symbol(x) || is_gap(x);
    else return sym_grp[CANONICAL].has_symbol(x);
}

bool Alphabet::is_gap(const char& x) const {
    return sym_grp[GAP].has_symbol(x);
}

bool Alphabet::is_degenerate(const char& x) const {
    return sym_grp[DEGENERATE].has_symbol(x);
}

bool Alphabet::is_unknown(const char& x) const {
    return sym_grp[UNKNOWN].has_symbol(x);
}

bool Alphabet::is_none(const char& x) const {
    return sym_grp[NONE].has_symbol(x);
}

bool Alphabet::is_missing(const char& x) const {
    return sym_grp[MISSING].has_symbol(x);
}

bool Alphabet::is_valid(const char& x) const {
    for (size_t g_idx = 0; g_idx < NUM_GROUP; ++g_idx) {
        if (sym_grp[g_idx].has_symbol(x)) return true;
    }
    return false;
}

void Alphabet::update_sym_idx() {
    std::string str = get_valid_symbol();
    for (size_t i = 0; i < str.size(); ++i) {
        this->sym_idx[str[i]] = i;
    }
}

std::string Alphabet::get_valid_symbol() const {
    std::string s;
    for (size_t g_idx = 0; g_idx < NUM_GROUP; ++g_idx) s += sym_grp[g_idx].get_member();
    return s;
}

std::string Alphabet::get_degeneracy(const char& x, float* w) const {
    std::string s;
    if (is_canonical(x)) s += x;
    else if (is_degenerate(x)) {
        std::pair<DegeneracyMap::const_iterator, DegeneracyMap::const_iterator> p = degeneracy.equal_range(x);
        for (DegeneracyMap::const_iterator pos = p.first; pos != p.second; ++pos) s += (*pos).second;
    }
    if (w != NULL) {
        if (s.empty()) *w = 0.;
        else *w = 1. / (float) s.size();
    }
    return s;
}

VectorXf Alphabet::get_count(const char& x) const {
    VectorXf v = VectorXf::Zero(get_canonical_size());
    if (is_canonical(x)) {
        v(get_idx(x)) = 1;
    } else if (is_degenerate(x)) {
        std::string s = get_degeneracy(x);
        for (std::string::iterator pos = s.begin(); pos != s.end(); ++pos)
            v(get_idx(*pos)) = 1;
        v /= v.sum();
    }
    return v;
}

/**
   @class AminoAcid
 */

AminoAcid::AminoAcid(const bool& nocase, const bool& gapres) 
: Alphabet("ACDEFGHIKLMNPQRSTVWY", "-", "BJZOU", "X", "*", "~", nocase, gapres) {
    init();
}

AminoAcid::AminoAcid(const char* gap, const bool& nocase, const bool& gapres) 
: Alphabet("ACDEFGHIKLMNPQRSTVWY", gap, "BJZOU", "X", "*", "~", nocase, gapres) {
    init();
}

void AminoAcid::init() {
    set_degeneracy('B', 'N');
    set_degeneracy('B', 'D');
    set_degeneracy('J', 'I');
    set_degeneracy('J', 'L');
    set_degeneracy('Z', 'Q');
    set_degeneracy('Z', 'E');
    set_degeneracy('O', 'K');
    set_degeneracy('U', 'C');
}
