#include "seqio.h"

const std::string TRIM_SPACE = " \t\n\v";

inline std::string rtrim(std::string& s, const std::string& drop = TRIM_SPACE) {
    return s.erase(s.find_last_not_of(drop) + 1);
}

inline std::string ltrim(std::string& s, const std::string& drop = TRIM_SPACE) {
    return s.erase(0, s.find_first_not_of(drop));
}

inline std::string trim(std::string& s, const std::string& drop = TRIM_SPACE) {
    std::string r = rtrim(s, drop);
    return ltrim(r, drop);
}  

FastaParser::FastaParser(std::istream& is) : is(is) { 
    init_pos = is.tellg();
    init();
}

SeqRecord FastaParser::next() {
    SeqRecord r;
    r.id = buffer.substr(1, buffer.find(' ') - 1);
    r.desc = buffer.substr(1);
    std::string s;
    while (getline(is, s)) {
        s = trim(s);
        if (s[0] == '>') {
            set_buffer(s);
            return r;
        }
        r.seq += s;
    }
    buffer.clear();
    _has_next = false;
    return r;
}

void FastaParser::init() {
    std::string s;
    is.seekg(init_pos);
    while (getline(is, s)) {
        s = trim(s);
        if (s.empty()) continue;
        if (s[0] == '>') {
            set_buffer(s);
            return;
        } else {
            // something wrong
            // raise exception 
        }
    }
}
