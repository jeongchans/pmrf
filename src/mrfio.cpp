#include "mrfio.h"

using std::endl;
using std::fixed;
using std::setprecision;

EdgeIndexVector EdgeIndexImporter::import(std::istream& is) {
    size_t idx1, idx2;
    EdgeIndexVector eidxs;
    while (is >> idx1) {
        is >> idx2;
        idx1--;
        idx2--;
        eidxs.push_back(EdgeIndex(idx1, idx2));
    }
    return eidxs;
}

void MRFExporter::export_model(const MRF& model, ostream& os) {
    export_header(model, os);
    export_body(model, os);
}

void MRFExporter::export_header(const MRF& model, ostream& os) {
    os << "LENG\t" << model.get_length() << endl;
    os << "SEQ" << endl;
    export_seq(model.get_seq(), 100, os);
}

void MRFExporter::export_body(const MRF& model, ostream& os) {
    string symbols = model.get_var_symbol();
    string seq = model.get_seq();
    size_t n = model.get_length();
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    os << "# PROFILE" << endl;
    os << "RES\t";
    export_node_symbol(symbols, os);
    os << endl;
    for (size_t i = 0; i < n; ++i) {
        os << seq[i] << " " << i + 1 << "\t";
        export_psfm(model.get_psfm(i), os);
        os << endl;
    }
    os << "# NODE" << endl;
    os << "RES\t";
    export_node_symbol(symbols, os);
    os << endl;
    for (size_t i = 0; i < n; ++i) {
        os << seq[i] << " " << i + 1 << "\t";
        export_node_weight(model.get_node(i).get_weight(), os);
        os << endl;
    }
    os << "# EDGE" << endl;
    os << "RES1\tRES2\t";
    export_edge_symbol(symbols, os);
    os << endl;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        size_t i = pos->idx1;
        size_t j = pos->idx2;
        os << seq[i] << " " << i + 1 << "\t";
        os << seq[j] << " " << j + 1 << "\t";
        export_edge_weight(model.get_edge(i, j).get_weight(), os);
        os << endl;
    }
    os << "//" << endl;
}

void MRFExporter::export_seq(const string& seq, const size_t& width, ostream& os) {
    size_t n = seq.size();
    size_t i = 0;
    do {
        os << seq.substr(i, width) << endl;
        i += width;
    } while (i < n);
}

void MRFExporter::export_node_symbol(const string& sym, ostream& os) {
    export_elem(sym[0], os);
    for (size_t i = 1; i < sym.size(); ++i) export_elem(sym[i], os, true);
}

void MRFExporter::export_psfm(const Float1dArray& w, ostream& os) {
    export_elem(w(0), os);
    for (int i = 1; i < w.size(); ++i) export_elem(w(i), os, true);
}

void MRFExporter::export_node_weight(const Float1dArray& w, ostream& os) {
    export_elem(w(0), os);
    for (int i = 1; i < w.size(); ++i) export_elem(w(i), os, true);
}

void MRFExporter::export_edge_symbol(const string& sym, ostream& os) {
    string s = "";
    s = sym[0];
    s += sym[0];
    export_elem(s, os);
    for (size_t i = 0; i < sym.size(); ++i) {
        for (size_t j = 0; j < sym.size(); ++j) {
            if (i == 0 && j == 0) continue;
            s = sym[i];
            s += sym[j];
            export_elem(s, os, true);
        }
    }
}

void MRFExporter::export_edge_weight(const Float2dArray& w, ostream& os) {
    export_elem(w(0, 0), os);
    for (int i = 0; i < w.rows(); ++i) {
        for (int j = 0; j < w.cols(); ++j) {
            if (i == 0 && j == 0) continue;
            export_elem(w(i, j), os, true);
        }
    }
}

void MRFExporter::export_elem(const double& x, ostream& os, bool pre_sep) {
    if (pre_sep) os << "\t";
    if (x == 0) os << "*";
    else os << setprecision(4) << x;
}

void MRFExporter::export_elem(const char& x, ostream& os, bool pre_sep) {
    if (pre_sep) os << "\t";
    os << x;
}

void MRFExporter::export_elem(const string& x, ostream& os, bool pre_sep) {
    if (pre_sep) os << "\t";
    os << x;
}

MRF MRFImporter::import_model(istream& is) {
    //TODO
}
