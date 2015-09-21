#include "trace.h"

#include <stdexcept>

#include "seq/seqio.h"

/**
   @class Trace
 */

//const char Trace::STATE_SYMBOL[5] = {'M', 'I', 'O', 'E', 'U'};

const char Trace::STATE_SYMBOL[3] = {'M', 'I', 'D'};

inline Trace::StateType Trace::char_to_state(const char& x) {
    for (size_t i = 0; i < NUM_STATE_TYPE; ++i)
        if (x == STATE_SYMBOL[i]) return (StateType)i;
    return UNDEFINED;
}

Trace::Trace(const string& ststr, const string& seq) {
    init(ststr, seq, "", "");
}

//Trace::Trace(const string& ststr, const string& seq, const string& id, const string& desc) : id(id), desc(desc) {
//    size_t i = 0;
//    size_t j = 0;
//    for (string::const_iterator pos = ststr.begin(); pos != ststr.end(); ++pos) {
//        char aa;
//        StateType st = char_to_state(*pos);
//        if (st == MATCH || st == INSERT) aa = seq[i++];
////        else if (st == GAPOPEN) aa = '=';
////        else if (st == GAPEXT) aa = '-';
////        else if (st == UNALIGNED) aa = '^';
//        else throw std::range_error(ststr);
//        this->st.push_back(State(st, aa));
//        if (st != INSERT) this->aidx.push_back(j);
//        ++j;
//    }
//}

Trace::Trace(const string& ststr, const string& seq, const string& id, const string& desc) {
    init(ststr, seq, id, desc);
}

void Trace::init(const string& ststr, const string& seq, const string& id, const string& desc) {
    this->id = id;
    this->desc = desc;
    size_t i = 0;
    size_t j = 0;
    for (string::const_iterator pos = ststr.begin(); pos != ststr.end(); ++pos) {
        char aa;
        StateType st = char_to_state(*pos);
        if (st == MATCH || st == INSERT) aa = seq[i++];
        else if (st == DELETE) aa = '-';
        else throw std::range_error(ststr);
        this->st.push_back(State(st, aa));
        if (st != INSERT) this->aidx.push_back(j);
        ++j;
    }
}

string Trace::get_seq() const {
    string seq;
    for (vector<State>::const_iterator pos = st.begin(); pos != st.end(); ++pos)
        if (pos->st == MATCH || pos->st == INSERT) seq += pos->aa;
    return seq;
}

string Trace::get_matched_aseq() const {
    string seq;
    for (vector<size_t>::const_iterator pos = aidx.begin(); pos != aidx.end(); ++pos)
        seq += st[*pos].aa;
    return seq;
}

bool Trace::operator==(const Trace& rhs) const {
    size_t n = st.size();
    if (n != rhs.st.size()) return false;
    for (size_t i = 0; i < n; ++i)
        if (st[i] != rhs.st[i]) return false;
    return true;
}

/**
   @class TraceVector
 */

vector<string> TraceVector::get_matched_aseq_vec() const {
    vector<string> r;
    for (const_iterator pos = begin(); pos != end(); ++pos)
        r.push_back(pos->get_matched_aseq());
    return r;
}

/**
   @class TraceImporter
 */

std::pair<size_t, TraceVector> TraceImporter::import(std::istream& is, const MSAFormat& fmt) {
    TraceVector traces;
    switch(fmt) {
    case AFASTA:
        afa_to_trace(is, traces);
        break;
    case A3M:
        a3m_to_trace(is, traces);
        break;
    }
    size_t length = 0;
    if (! traces.empty()) length = traces[0].get_length();
    return make_pair(length, traces);
}

//inline string reformat_a3m_seq(string aseq) {
//    for (string::iterator pos = aseq.begin(); pos != aseq.end(); ++pos) {
//        if (*pos == '-') *pos = '^';
//        else break;
//    }
//    for (string::reverse_iterator pos = aseq.rbegin(); pos != aseq.rend(); ++pos) {
//        if (*pos == '-') *pos = '^';
//        else break;
//    }
//    bool gap_ahead = false;
//    for (string::iterator pos = aseq.begin(); pos != aseq.end(); ++pos) {
//        if (*pos == '-') {
//            if (!gap_ahead) *pos = '=';
//            gap_ahead = true;
//        } else {
//            gap_ahead = false;
//        }
//    }
//    return aseq;
//}

inline string reformat_a3m_seq(string aseq) {
    return aseq;
}

string reformat_afa_seq(string aseq, const string& ref_aseq) {
    string a3m_seq;
    size_t n = ref_aseq.size();
    for (size_t i = 0; i < n; ++i) {
        if (ref_aseq[i] != '-') a3m_seq += toupper(aseq[i]);
        else if (aseq[i] != '-') a3m_seq += tolower(aseq[i]);
    }
    return reformat_a3m_seq(a3m_seq);
}

//void TraceImporter::parse_state(const std::string& aseq, std::string& st, std::string& seq) {
//    st.empty();
//    seq.empty();
//    for (string::const_iterator pos = aseq.begin(); pos != aseq.end(); ++pos) {
//        if (*pos == '=') st += 'O';
//        else if (*pos == '-') st += 'E';
//        else if (*pos == '^') st += 'U';
//        else if (abc.is_valid(toupper(*pos))) {
//            if (islower(*pos)) st += 'I';
//            else st += 'M';
//            seq += toupper(*pos);
//        } else {
//            std::cerr << "Undefined symbol: " << *pos << std::endl;
//            throw std::range_error(aseq);
//        }
//    }
//}

void TraceImporter::parse_state(const std::string& aseq, std::string& st, std::string& seq) {
    st.empty();
    seq.empty();
    for (string::const_iterator pos = aseq.begin(); pos != aseq.end(); ++pos) {
        if (*pos == '-') st += 'D';
        else if (abc.is_valid(toupper(*pos))) {
            if (islower(*pos)) st += 'I';
            else st += 'M';
            seq += toupper(*pos);
        } else {
            std::cerr << "Undefined symbol: " << *pos << std::endl;
            throw std::range_error(aseq);
        }
    }
}

bool TraceImporter::is_valid_state(const std::string& st) const {
    return st.find_first_of('M') != string::npos;
}

void TraceImporter::a3m_to_trace(std::istream& is, TraceVector& traces) {
    FastaParser parser(is);
    while (parser.has_next()) {
        SeqRecord r = parser.next();
        r.seq = reformat_a3m_seq(r.seq);
        string st, seq;
        parse_state(r.seq, st, seq);
        if (is_valid_state(st)) traces.push_back(Trace(st, seq, r.id, r.desc));
    }
}

void TraceImporter::afa_to_trace(std::istream& is, TraceVector& traces) {
    string ref_aseq;
    FastaParser parser(is);
    if (parser.has_next()) {
        SeqRecord r = parser.next();
        ref_aseq = r.seq;
    }
    parser.init();
    while (parser.has_next()) {
        SeqRecord r = parser.next();
        r.seq = reformat_afa_seq(r.seq, ref_aseq);
        string st, seq;
        parse_state(r.seq, st, seq);
        if (is_valid_state(st)) traces.push_back(Trace(st, seq, r.id, r.desc));
    }
}
