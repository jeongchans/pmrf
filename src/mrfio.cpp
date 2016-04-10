#include "mrfio.h"

#include <iomanip>

#include "command.h"

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
    const string sep = "\t";
    Eigen::IOFormat fmt(4, Eigen::DontAlignCols, sep, sep);
    size_t n = model.get_length();
    string seq = model.get_seq();
    size_t k = model.get_num_var();
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    size_t n2 = edge_idxs.size();
    os.write((char*)&n, sizeof(size_t));
    os.write(seq.c_str(), seq.size() + 1);
    // Writing MRF architecture
    os.write((char*)&n2, sizeof(size_t));
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        os.write((char*)&(pos->idx1), sizeof(size_t));
        os.write((char*)&(pos->idx2), sizeof(size_t));
    }
    // Writing sequence profile
    os.write((char*)model.get_psfm().data(), n * k * sizeof(float));
    // Writing node weight
    for (size_t i = 0; i < n; ++i)
        os.write((char*)model.get_node(i).get_weight().data(), k * sizeof(float));
    // Writing edge weight
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos)
        os.write((char*)model.get_edge(pos->idx1, pos->idx2).get_weight().data(), k * k * sizeof(float));
}

MRF MRFImporter::import_model(istream& is, const Alphabet& abc) {
    size_t n, n2;
    string seq;
    is.read((char*)&n, sizeof(size_t));
    std::getline(is, seq, '\0');
    // Reading MRF architecture
    is.read((char*)&n2, sizeof(size_t));
    EdgeIndexVector edge_idxs;
    size_t idx1, idx2;
    for (size_t i = 0; i < n2; ++i) {
        is.read((char*)&idx1, sizeof(size_t));
        is.read((char*)&idx2, sizeof(size_t));
        edge_idxs.push_back(EdgeIndex(idx1, idx2));
    }
    MRF model(seq, abc, &edge_idxs);
    size_t k = model.get_num_var();
    // Reading profile
    MatrixXf m(n, k);
    is.read((char*)m.data(), n * k * sizeof(float));
    model.set_psfm(m);
    // Reading node weight
    VectorXf v(k);
    for (size_t i = 0; i < n; ++i) {
        is.read((char*)v.data(), k * sizeof(float));
        model.get_node(i).set_weight(v);
    }
    // Reading edge weight
    m.resize(k, k);
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        is.read((char*)m.data(), k * k * sizeof(float));
        model.get_edge(pos->idx1, pos->idx2).set_weight(m);
    }
    return model;
}

string format_seq(const string& seq, const size_t& width) {
    string ret = seq.substr(0, width);
    for (size_t i = width; i < seq.size(); i += width) ret += "\n" + seq.substr(i, width);
    return ret;
}

string format_symbol(const string& symbol, const string& sep) {
    string ret;
    ret += symbol[0];
    for (size_t i = 1; i < symbol.size(); ++i) ret += sep + symbol[i];
    return ret;
}

string format_paired_symbol(const string& symbol, const string& sep) {
    string ret;
    ret += symbol[0];
    ret += symbol[0];
    for (size_t i = 0; i < symbol.size(); ++i)
        for (size_t j = 0; j < symbol.size(); ++j)
            if (i != 0 || j != 0) ret += sep + symbol[i] + symbol[j];
    return ret;
}

ostream& operator<<(ostream& os, const MRF& model) {
    const string sep = "\t";
    Eigen::IOFormat fmt(4, Eigen::DontAlignCols, sep, sep);
    size_t n = model.get_length();
    string seq = model.get_seq();
    string symbol = model.get_var_symbol();
    EdgeIndexVector edge_idxs = model.get_edge_idxs();
    os << "VER" << sep << VERSION << endl
       << "LENGTH" << sep << n << endl
       << "SEQ" << endl
       << format_seq(seq, 100) << endl;
    os << "# PROFILE" << endl
       << "RES" << sep << format_symbol(symbol, sep) << endl;
    for (size_t i = 0; i < n; ++i)
        os << seq[i] << " " << i + 1 << sep << model.get_psfm(i).format(fmt) << endl;
    os << "# NODE" << endl
       << "RES" << sep << format_symbol(symbol, sep) << endl;
    for (size_t i = 0; i < n; ++i)
        os << seq[i] << " " << i + 1 << sep << model.get_node(i).get_weight().format(fmt) << endl;
    os << "# EDGE" << endl
       << "RES1" << sep << "RES2" << sep << format_paired_symbol(symbol, sep) << endl;
    for (EdgeIndexVector::const_iterator pos = edge_idxs.begin(); pos != edge_idxs.end(); ++pos) {
        size_t i = pos->idx1;
        size_t j = pos->idx2;
        os << seq[i] << " " << i + 1 << sep
           << seq[j] << " " << j + 1 << sep
           << model.get_edge(i, j).get_weight().format(fmt) << endl;
    }
    os << "//" << endl;
    return os;
}
