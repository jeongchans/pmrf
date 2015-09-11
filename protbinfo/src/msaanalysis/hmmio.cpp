#include "hmmio.h"

//#include <iomanip>
//#include <cmath>
//
//using std::setw;
//using std::left;
using std::endl;
using std::fixed;
using std::setprecision;

void HHsearchHMMExporter::export_model(const ProfileHMM& model, std::ostream& os) {
    export_header(model, os);
    os << "#" << endl;
    export_body(model, os);
}

void HHsearchHMMExporter::export_header(const ProfileHMM& model, std::ostream& os) {
    string desc = model.get_desc();
    if (!desc.empty()) os << "NAME\t" << desc << endl;
    os << "LENG\t" << model.get_length() << endl
       << "NEFF\t" << model.get_eff_num() << endl;
    os << "SEQ" << endl;
    os << ">" << desc << endl;
    export_seq(model.get_seq(), 100, os);
}

void HHsearchHMMExporter::export_seq(const std::string& seq, const size_t& width, std::ostream& os) {
    size_t n = seq.size();
    size_t i = 0;
    do {
        os << seq.substr(i, width) << endl;
        i += width;
    } while (i < n);
}

void HHsearchHMMExporter::export_body(const ProfileHMM& model, std::ostream& os) {
    os << "NULL\t";
    export_emit(model.get_null_emit(), os);
    os << endl;
    os << "HMM\t";
    string symbols = model.get_emit_symbol();
    export_emit_symbol(symbols, os);
    os << endl
       << "\tM->M\tM->I\tM->D\tI->M\tI->I\tD->M\tD->D\tNeff\tNeff_I\tNeff_D" << endl;
    os << "\t";
    export_transit(MATCH, model.get_null_transit(MATCH), os);
    export_transit(INSERT, model.get_null_transit(INSERT), os);
    export_transit(DELETE, model.get_null_transit(DELETE), os);
    os << "*\t*\t*" << endl;
    size_t n = model.get_length();
    std::string seq = model.get_seq();
    for (size_t i = 0; i < n; ++i) {
        const HMMMatchState& m_state = model.get_match(i);
        const HMMDeleteState& d_state = model.get_delete(i);
        const HMMInsertState& i_state = model.get_insert(i);
        size_t idx = i + 1;
        os << seq[i] << " " << idx << "\t";
        export_emit(m_state.get_emit(), os);
        os << idx << endl
           << "\t";
        export_transit(MATCH, m_state.get_transit(), os);
        export_transit(INSERT, i_state.get_transit(), os);
        export_transit(DELETE, d_state.get_transit(), os);
        export_eff_num(m_state.get_eff_num(), os);
        export_eff_num(i_state.get_eff_num(), os);
        export_eff_num(d_state.get_eff_num(), os);
        os << endl << endl;
    }
    os << "//" << endl;
}

void HHsearchHMMExporter::export_emit_symbol(const std::string& sym, std::ostream& os) {
    for (size_t i = 0; i < sym.size(); ++i) export_elem(sym[i], os);
}

void HHsearchHMMExporter::export_emit(const Float1dArray& v, std::ostream& os) {
    for (int i = 0; i < v.size(); ++i) export_elem(v(i), os);
}

void HHsearchHMMExporter::export_transit(const StateType& type, const Float1dArray& v, std::ostream& os) {
    switch (type) {
    case MATCH:
        export_elem(v(MATCH), os);
        export_elem(v(INSERT), os);
        export_elem(v(DELETE), os);
        break;
    case DELETE:
        export_elem(v(MATCH), os);
        export_elem(v(DELETE), os);
        break;
    case INSERT:
        export_elem(v(MATCH), os);
        export_elem(v(INSERT), os);
        break;
    default:
        throw;
    }
}

void HHsearchHMMExporter::export_eff_num(const double& v, std::ostream& os) {
    export_elem(1. / (v + 1.), os);
}

inline void HHsearchHMMExporter::export_elem(const double& x, std::ostream& os) {
    if (x == 1) os << 0;    // to avoid -0 problem
    else if (x > 0) os << fixed << setprecision(0) << (-1000. * log2(x));
    else os << "*";
    os << "\t";
}

inline void HHsearchHMMExporter::export_elem(const char& x, std::ostream& os) {
    os << x << "\t";
}

ProfileHMM HHsearchHMMImporter::import_model(std::istream& is) {
    ProfileHMM model = import_header(is);
    import_body(model, is);
    return model;
}

ProfileHMM HHsearchHMMImporter::import_header(std::istream& is) {
    string dummy, s;
    size_t length;
    string seq;
    string desc;
    double eff_num;
    while (is) {
        is >> s;
        if (s == "#") break;
        else if (s == "NAME") {
            getline(is, desc);
            desc = desc.substr(1);
        }
        else if (s == "LENG") {
            is >> length;
            getline(is, dummy);
        }
        else if (s == "NEFF") {
            is >> eff_num;
            getline(is, dummy);
        }
        else if (s == "SEQ") {
            getline(is, dummy);
            getline(is, dummy);
            seq = import_seq(length, is);
        }
    }
    ProfileHMM model(length, abc);
    model.set_seq(seq);
    model.set_desc(desc);
    model.set_eff_num(eff_num);
    return model;
}

std::string HHsearchHMMImporter::import_seq(const size_t& length, std::istream& is) {
    string seq, s;
    while (is) {
        getline(is, s);
        seq += s;
        if (seq.size() >= length) break;
    }
    return seq;
}

void HHsearchHMMImporter::import_body(ProfileHMM& model, std::istream& is) {
    std::string s, dummy;
    size_t length = model.get_length();
    is >> s;
    int sym_num = abc.get_canonical_size();
    model.set_null_emit(import_emit(is, sym_num));
    is >> s;
    string emit_symbol = import_emit_symbol(is);
    assert(abc.get_canonical() == emit_symbol);
    getline(is, s);
    model.set_null_transit(MATCH, import_transit(is, MATCH));
    model.set_null_transit(INSERT, import_transit(is, INSERT));
    model.set_null_transit(DELETE, import_transit(is, DELETE));
    getline(is, s);
    for (size_t i = 0; i < length; ++i) {
        is >> s;
        is >> dummy;
        model.get_match(i).set_emit(import_emit(is, sym_num));
        is >> dummy;
        model.get_match(i).set_transit(import_transit(is, MATCH));
        model.get_insert(i).set_transit(import_transit(is, INSERT));
        model.get_delete(i).set_transit(import_transit(is, DELETE));
        model.get_match(i).set_eff_num(import_eff_num(is));
        model.get_insert(i).set_eff_num(import_eff_num(is));
        model.get_delete(i).set_eff_num(import_eff_num(is));
    }
    while (is) {
        is >> s;
        if (s == "//") break;
    }
}

string HHsearchHMMImporter::import_emit_symbol(std::istream& is) {
    string dummy;
    string s = "";
    getline(is, dummy);
    for (std::string::iterator pos = dummy.begin(); pos != dummy.end(); ++pos)
        if (*pos != '\t' and *pos != ' ') s += *pos;
    return s;
}

Float1dArray HHsearchHMMImporter::import_emit(std::istream& is, const int& size) {
    Float1dArray v(size);
    for (int i = 0; i < size; ++i) v(i) = import_elem(is);
    return v;
}

Float1dArray HHsearchHMMImporter::import_transit(std::istream& is, const StateType& type) {
    Float1dArray v(NUM_STATE_TYPE);
    v = 0;
    switch (type) {
    case MATCH:
        v(MATCH) = import_elem(is);
        v(INSERT) = import_elem(is);
        v(DELETE) = import_elem(is);
        break;
    case DELETE:
        v(MATCH) = import_elem(is);
        v(DELETE) = import_elem(is);
        break;
    case INSERT:
        v(MATCH) = import_elem(is);
        v(INSERT) = import_elem(is);
        break;
    default:
        throw;
    }
    return v;
}

double HHsearchHMMImporter::import_eff_num(std::istream& is) {
    double d;
    is >> d;
    return 1. / pow(2., d / -1000.) - 1.;
}

inline double HHsearchHMMImporter::import_elem(std::istream& is) {
    string s;
    is >> s;
    if (s == "*") return 0.;
    else {
        double d;
        std::istringstream(s) >> d;
        return pow(2., d / -1000.);
    }
}

//HCIModel HCIModelImporter::import_model(std::istream& is) {
//    std::string dummy, s;
//    size_t length;
//    string desc, ss_dssp, ss_pred, ss_conf;
//    double eff_num;
//    while (is) {
//        is >> s;
//        if (s == "#") break;
//        else if (s == "NAME") {
//            getline(is, desc);
//            desc = desc.substr(1);
//        }
//        else if (s == "LENG") is >> length;
//        else if (s == "NEFF") is >> eff_num;
//        else if (s == "SS_DSSP") is >> ss_dssp;
//        else if (s == "SS_PRED") is >> ss_pred;
//        else if (s == "SS_CONF") is >> ss_conf;
//    }
//    if (is.eof()) throw NoEntryException();
//    getline(is, s);
//    HCIModel model(length, abc);
//    model.set_desc(desc);
//    model.set_eff_num(eff_num);
//    model.set_sec_struct(ss_dssp, ss_pred, ss_conf);
//    is >> s;
//    int sym_num = abc.get_canonical_size();
//    model.set_null_emit(import_emit(is, sym_num));
//    model.set_null_transit(MATCH, import_transit(is, MATCH));
//    model.set_null_transit(INSERT, import_transit(is, INSERT));
//    model.set_null_transit(DELETE, import_transit(is, DELETE));
//    getline(is, s);
//    is >> s;
//    std::string emit_symbol = import_emit_symbol(is);
//    assert(abc.get_canonical() == emit_symbol);
//    getline(is, s);
//    for (size_t i = 0; i < emit_symbol.size(); ++i) getline(is, s);
//    std::string seq = "";
//    for (int i = 0; i < length; ++i) {
//        is >> s;
//        seq += s;
//        is >> dummy;
//        model.set_emit(MATCH, i, import_emit(is, sym_num));
//        is >> dummy;
//        model.set_transit(MATCH, i, import_transit(is, MATCH));
//        model.set_transit(INSERT, i, import_transit(is, INSERT));
//        model.set_transit(DELETE, i, import_transit(is, DELETE));
//        model.set_eff_num(MATCH, i, import_eff_num(is));
//        model.set_eff_num(INSERT, i, import_eff_num(is));
//        model.set_eff_num(DELETE, i, import_eff_num(is));
//        model.set_property(MATCH, i, "cm_wtavg", import_cm_avg(is, sym_num));
//    }
//    model.set_seq(seq);
//
//    while (is) {
//        is >> s;
//        if (s == "//") break;
//        else if (s == "HCISVM_START") {
//            string hcisvm_content = extract_hcisvm_content(is);
//            std::istringstream iss(hcisvm_content);
//            HCISVMModelImporter hcisvm_importer;
//            HCISVMModel *hcisvm = new HCISVMModel;
//            hcisvm_importer.import_model(*hcisvm, iss);
//            model.set_hcisvm_model(hcisvm);
//        }
//    }
//
//    return model;
//}
//
//string HCIModelImporter::extract_hcisvm_content(std::istream& is) {
//    string dummy;
//    string s = "";
//    while (true) {
//        getline(is, dummy);
//        if (dummy.find("HCISVM_END") == 0) break;
//        s += dummy + "\n";
//    }
//    return s;
//}
//
//Float2dArray HCIModelImporter::import_cm_avg(std::istream& is, const int& size) {
//    Float2dArray x(size, size);
//    for (int i = 0; i < size; ++i) {
//        for (int j = 0; j < size; ++j) x(i, j) = -log(1. / import_elem(is) - 1.);
//    }
//    return x;
//}
//
//std::string HCIModelImporter::import_emit_symbol(std::istream& is) {
//    std::string dummy;
//    std::string s = "";
//    getline(is, dummy);
//    for (std::string::iterator pos = dummy.begin(); pos != dummy.end(); ++pos)
//        if (*pos != '\t') s += *pos;
//    return s;
//}
//
//Float1dArray HCIModelImporter::import_emit(std::istream& is, const int& size) {
//    Float1dArray v(size);
//    for (int i = 0; i < size; ++i) v(i) = import_elem(is);
//    return v;
//}
//
//Float1dArray HCIModelImporter::import_transit(std::istream& is, const StateType& type) {
//    Float1dArray v(NUM_STATE_TYPE);
//    v = 0;
//    switch (type) {
//    case MATCH:
//        v(MATCH) = import_elem(is);
//        v(INSERT) = import_elem(is);
//        v(DELETE) = import_elem(is);
//        break;
//    case DELETE:
//        v(MATCH) = import_elem(is);
//        v(DELETE) = import_elem(is);
//        break;
//    case INSERT:
//        v(MATCH) = import_elem(is);
//        v(INSERT) = import_elem(is);
//        break;
//    }
//    return v;
//}
//
//double HCIModelImporter::import_eff_num(std::istream& is) {
//    double d;
//    is >> d;
//    return 1. / pow(2., d / -1000.) - 1.;
//}
//
//inline double HCIModelImporter::import_elem(std::istream& is) {
//    std::string s;
//    is >> s;
//    if (s == "*") return 0.;
//    return pow(2., atof(s.c_str()) / -1000.);
//}
