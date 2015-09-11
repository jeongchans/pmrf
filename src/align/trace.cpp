#include "trace.h"

#include <algorithm>

#include "seq/seqio.h"

inline StateType char_to_state_type(const char& x) {
    if (x == 'M') return MATCH;
    else if (x == 'D') return DELETE;
    else if (x == 'I') return INSERT;
    else return UNDEFINED;
}

Trace::Trace() {
}

Trace::Trace(const char* st, const char* seq) {
    init_trace(string(st), string(seq));
}

Trace::Trace(const string& st, const string& seq) {
    init_trace(st, seq);
}

Trace::Trace(const string& st, const string& seq, const string& id, const string& desc) : id(id), desc(desc) {
    init_trace(st, seq);
}

bool Trace::operator==(const Trace& rhs) const {
    for (size_t pos = 0; pos < length; ++pos) {
        for (size_t type = 0; type < NUM_STATE_TYPE; ++type) {
            if (emit[type][pos] != rhs.emit[type][pos]) return false;
        }
    }
    return true;
}

string Trace::get_seq() const {
    string seq = "";
    for (size_t i = 0; i < length; ++i)
        seq += emit[MATCH][i] + emit[INSERT][i];
    return seq;
}

bool Trace::is_passing(const StateType& type, const size_t& pos) const {
    if (emit[type][pos].empty()) return false;
    return true;
}

std::string Trace::get_emit(const StateType& type, const size_t& pos) const {
    return emit[type][pos];
}

std::string Trace::get_MD_seq() const {
    std::string s;
    for (size_t i = 0; i < length; ++i) {
        s += emit[MATCH][i] + emit[DELETE][i];
    }
    return s;
}

bool Trace::has_terminal_gap(const size_t& pos) const {
    if (pos < first_not_delete || pos > last_not_delete) return true;
    return false;
}

int Trace::get_visit(const StateType& type, const size_t& pos) const {
    return (int) emit[type][pos].size();
}

FreqVec Trace::get_transit_count(const StateType& type, const size_t& pos) const {
    FreqVec fv(NUM_STATE_TYPE);
    fv = 0;
    if (is_passing(type, pos)) {
        if (pos == length - 1) {
            if (type == MATCH && get_visit(INSERT, pos) > 0) fv(INSERT) = 1;
        }
        else if (type == MATCH) {
            if (get_visit(INSERT, pos) > 0) fv(INSERT) = 1;
            else if (get_visit(MATCH, pos + 1) > 0) fv(MATCH) = 1;
            else fv(DELETE) = 1;
        }
        else if (type == DELETE) {
            if (get_visit(MATCH, pos + 1) > 0) fv(MATCH) = 1;
            else fv(DELETE) = 1;
        }
        else if (type == INSERT) {
            fv(INSERT) = get_visit(INSERT, pos) - 1;
            fv(MATCH) = 1;
        }
    }
    return fv;
}

void Trace::init_trace(const string& st, const string& seq) {
    length = count(st.begin(), st.end(), 'M') + count(st.begin(), st.end(), 'D');
    for (size_t i = 0; i < NUM_STATE_TYPE; ++i) {
        emit[i] = std::vector<std::string>(length, "");
    }
    int idx = -1;
    StateType cur_type;
    int seq_idx = -1;
    std::string::const_iterator pos = st.begin();
    while (true) {
        cur_type = char_to_state_type(*pos);
        if (cur_type == INSERT && idx < 0) {
            ++seq_idx;
            if (++pos == st.end()) break;
            continue;
        }
        if (cur_type == MATCH || cur_type == DELETE) ++idx;
        if (cur_type == MATCH || cur_type == INSERT) emit[cur_type][idx] += seq[++seq_idx];
        else emit[DELETE][idx] += "-";
        if (++pos == st.end()) break;
    }
    // identify terminal gap
    first_not_delete = st.find_first_not_of('D');
    last_not_delete = length - (st.size() - st.find_last_not_of('D'));
}

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

void TraceImporter::afa_to_trace(std::istream& is, TraceVector& traces) {
    vector<SeqRecord> records;
    FastaParser parser(is);
    while (parser.has_next()) {
        SeqRecord r = parser.next();
        if (r.id == "ss_dssp") traces.ss_dssp = r.seq;
        else if (r.id == "ss_pred") traces.ss_pred = r.seq;
        else if (r.id == "ss_conf") traces.ss_conf = r.seq;
        else records.push_back(r);
    }
    SeqRecord& ref = records[0];
    size_t n = ref.seq.size();
    for (vector<SeqRecord>::iterator pos = records.begin(); pos != records.end(); ++pos) {
        string st, seq;
        for (size_t i = 0; i < n; ++i) {
            char c = toupper(pos->seq[i]);
            char ref_c = toupper(ref.seq[i]);
            if (abc.is_gap(c)) {
                if (ref_c == '-') continue;
                st += 'D';
            } else if (abc.is_valid(c)) {
                if (ref_c == '-' && c == '-') continue;
                else if (ref_c == '-' && c != '-') st += 'I';
                else st += 'M';
                seq += c;
            } else {
                std::cerr << "Undefined symbol: " << pos->seq[i] << std::endl;
                // undefined symbol
                // raise exception
            }
        }
        Trace t(st, seq, pos->id, pos->desc);
        if (st.find_first_of('M') != string::npos && t.get_MD_seq().find_first_not_of('-') != string::npos) {
            traces.push_back(t);
        }
    }
}

void TraceImporter::a3m_to_trace(std::istream& is, TraceVector& traces) {
    FastaParser parser(is);
    while (parser.has_next()) {
        SeqRecord r = parser.next();
        if (r.id == "ss_dssp") {
            traces.ss_dssp = r.seq;
            continue;
        } else if (r.id == "ss_pred") {
            traces.ss_pred = r.seq;
            continue;
        } else if (r.id == "ss_conf") {
            traces.ss_conf = r.seq;
            continue;
        }
        std::string st, seq;
        for (std::string::iterator pos = r.seq.begin(); pos != r.seq.end(); ++pos) {
            if (abc.is_gap(*pos)) {
                st += 'D';
            } else if (abc.is_valid(toupper(*pos))) {
                if (islower(*pos)) st += 'I';
                else st += 'M';
                seq += toupper(*pos);
            } else {
                std::cerr << "Undefined symbol: " << *pos << std::endl;
                // undefined symbol
                // raise exception
            }
        }
        if (st.find_first_of('M') != string::npos) {
            traces.push_back(Trace(st, seq, r.id, r.desc));
        }
    }
}

TraceVector TraceVector::subset_passing(const StateType& type, const size_t& idx, const bool& omit_termi_gap) const {
    TraceVector trs;
    for (const_iterator pos = begin(); pos != end(); ++pos) {
        if (pos->is_passing(type, idx)) {
            if (omit_termi_gap && pos->has_terminal_gap(idx)) continue;
            trs.push_back(*pos);
        }
    }
    return trs;
}

vector<string> TraceVector::get_MD_seq_vec() const {
    vector<string> r;
    for (const_iterator pos = begin(); pos != end(); ++pos)
        r.push_back(pos->get_MD_seq());
    return r;
}
