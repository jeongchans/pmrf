#include "mrfio.h"

EdgeIndexVector EdgeIndexImporter::import(istream& is) {
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
    // Writing modeling log
    size_t fmtnum = model.get_fmtnum();
    os.write((char*)&fmtnum, sizeof(size_t));
    float neff = model.get_neff();
    os.write((char*)&neff, sizeof(float));
    //string ver = VERSION;
    // Writing sequence
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
    // Reading logging information
    size_t fmtnum;
    is.read((char*)&fmtnum, sizeof(size_t));
    float neff;
    is.read((char*)&neff, sizeof(float));
    // Reading sequence
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
    model.set_fmtnum(fmtnum);
    model.set_neff(neff);
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
